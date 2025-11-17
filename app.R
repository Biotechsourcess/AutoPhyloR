# --- Carga de Librerías ---
# Estas son las dependencias de la *aplicación*, no del paquete.
library(shiny)
library(shinyjs) # Para habilitar/deshabilitar pestañas
library(future)
library(promises)
library(ggtree)   # Para renderPlot
library(ggplot2)  # Para renderPlot y ggsave
library(Biostrings) # Para writeXStringSet
library(ape)        # Para write.tree

# Configura 'future' para usar procesos de R en segundo plano (multisesión)
future::plan(multisession)

# --- Definición de la Interfaz de Usuario (UI) ---
ui <- fluidPage(
  useShinyjs(), # Inicializar shinyjs
  titlePanel("AutoPhyloR: Pipeline Filogenético"),
  
  sidebarLayout(
    sidebarPanel(
      # Usar conditionalPanel para mostrar controles específicos del paso
      
      # --- PASO 1: ENTRADA ---
      conditionalPanel(
        condition = "input.wizardTabs == 'tab_1_input'",
        h4("Paso 1: Obtener Secuencias"),
        
        # Opción A: Búsqueda en NCBI (de su script)
        h5("Opción A: Búsqueda en NCBI"),
        textInput("searchTerm", "Término de Búsqueda (Gen/Proteína):", "COI"),
        textAreaInput("speciesVector", "Especies (una por línea):", "Homo sapiens\nPan troglodytes\nMus musculus\nGorilla gorilla"),
        selectInput("dbType", "Tipo de Base de Datos:", c("Proteína" = "protein", "Nucleótido" = "nuccore")),
        textInput("email", "Su Email (Requerido por NCBI):", ""),
        actionButton("run_ncbi_search", "1. Iniciar Búsqueda NCBI", class = "btn-primary"),
        
        hr(),
        
        # Opción B: Carga de FASTA (del documento)
        h5("Opción B: Cargar archivo FASTA"),
        fileInput("fasta_file_upload", "Cargar archivo FASTA", accept = c(".fasta", ".fa", ".fas")),
        selectInput("seqType_upload", "Tipo de Secuencia (si carga archivo):", c("DNA" = "DNA", "Proteína" = "AA")),
        actionButton("run_file_upload", "1. Cargar Secuencias", class = "btn-primary")
      ),
      
      # --- PASO 2: MSA ---
      conditionalPanel(
        condition = "input.wizardTabs == 'tab_2_msa'",
        h4("Paso 2: Alineamiento Múltiple"),
        p("Las secuencias se han cargado. Seleccione el método de MSA."),
        selectInput("msaMethod", "Método de MSA:", 
                    choices = c("ClustalOmega", "ClustalW", "Muscle"), #
                    selected = "ClustalOmega"),
        actionButton("run_msa", "2. Ejecutar MSA", class = "btn-primary")
      ),
      
      # --- PASO 3: INFERENCIA ---
      conditionalPanel(
        condition = "input.wizardTabs == 'tab_3_tree'",
        h4("Paso 3: Inferencia del Árbol"),
        p("El alineamiento está completo. Seleccione el método de inferencia."),
        selectInput("phyloMethod", "Método de Inferencia:",
                    choices = c("Neighbor-Joining (Distancia)" = "NJ", #
                                "Máxima Parsimonia" = "Parsimonia",         #
                                "Máxima Verosimilitud (Rápida)" = "ML"),    #
                    selected = "NJ"),
        actionButton("run_tree", "3. Generar Árbol", class = "btn-primary")
      ),
      
      # --- PASO 4: RESULTADOS ---
      conditionalPanel(
        condition = "input.wizardTabs == 'tab_4_results'",
        h4("Paso 4: Descargar Resultados"),
        p("Pipeline completado. Puede descargar los artefactos."),
        # Botones de descarga
        downloadButton("download_tree", "Descargar Árbol (Newick)"),
        downloadButton("download_msa", "Descargar Alineamiento (FASTA)"),
        downloadButton("download_plot", "Descargar Imagen del Árbol (PNG)")
      )
    ),
    
    mainPanel(
      # El "wizard" usa un tabsetPanel
      tabsetPanel(
        id = "wizardTabs",
        
        # Pestaña 1: Entrada
        tabPanel("1. Entrada de Secuencias", value = "tab_1_input",
                 h4("Esperando entrada de secuencias..."),
                 p("Complete los campos en el panel izquierdo (Opción A o B) y presione el botón correspondiente para comenzar.")
        ),
        
        # Pestaña 2: MSA
        tabPanel("2. Alineamiento (MSA)", value = "tab_2_msa",
                 h4("Alineamiento Múltiple de Secuencias"),
                 p("Presione 'Ejecutar MSA' para iniciar el alineamiento. Esto puede tomar varios minutos."),
                 verbatimTextOutput("msa_print_output") #
        ),
        
        # Pestaña 3: Árbol
        tabPanel("3. Inferencia de Árbol", value = "tab_3_tree",
                 h4("Inferencia de Árbol Filogenético"),
                 p("Presione 'Generar Árbol' para iniciar la inferencia. Esto puede ser rápido (NJ) o lento (ML/Parsimonia).")
        ),
        
        # Pestaña 4: Resultados
        tabPanel("4. Resultados y Visualización", value = "tab_4_results",
                 h4("Visualización del Árbol"),
                 plotOutput("tree_plot_output", height = "600px"), #
                 hr(),
                 h4("Alineamiento (Vista previa)"),
                 verbatimTextOutput("msa_final_summary")
        )
      )
    )
  )
)

# --- Definición del Servidor ---
server <- function(input, output, session) {
  
  # --- Gestión del Estado del Pipeline ---
  pipeline_state <- reactiveValues(
    raw_sequences = NULL,
    sequence_type = "DNA",
    msa_alignment = NULL,
    phylo_tree = NULL
  )
  
  # --- Deshabilitar pestañas futuras al inicio ---
  observe({
    shinyjs::disable(selector = "#wizardTabs li a[data-value='tab_2_msa']")
    shinyjs::disable(selector = "#wizardTabs li a[data-value='tab_3_tree']")
    shinyjs::disable(selector = "#wizardTabs li a[data-value='tab_4_results']")
  })
  
  
  # --- Lógica del Botón: Paso 1A (Búsqueda NCBI) ---
  observeEvent(input$run_ncbi_search, {
    
    # Validaciones
    if (is.null(input$email) || input$email == "") {
      showModal(modalDialog(title = "Error", "Se requiere un email válido para usar la API de NCBI."))
      return()
    }
    species_list <- unlist(strsplit(input$speciesVector, "\n"))
    species_list <- species_list[species_list != ""]
    if (length(species_list) < 3) {
      showModal(modalDialog(title = "Error", "Se requieren al menos 3 especies."))
      return()
    }
    
    showModal(modalDialog("Iniciando búsqueda en NCBI...", title = "Procesando", footer = NULL))
    
    # Copiar variables
    sTerm <- input$searchTerm
    sVec <- species_list
    dType <- input$dbType
    sEmail <- input$email
    
    # --- Ejecución Asíncrona ---
    future({
      AutoPhyloR::run_ncbi_search(searchTerm = sTerm,
                                  speciesVector = sVec,
                                  dbType = dType,
                                  email = sEmail)
      # SINTAXIS CORREGIDA A (function(...) { ... })
    }) %...>% (function(result_sequences) {
      # --- Éxito ---
      removeModal()
      pipeline_state$raw_sequences <- result_sequences
      pipeline_state$sequence_type <- if (dType == "protein") "AA" else "DNA"
      shinyjs::enable(selector = "#wizardTabs li a[data-value='tab_2_msa']")
      updateTabsetPanel(session, "wizardTabs", selected = "tab_2_msa")
      
    }) %...>% (function(error) {
      # --- Error (CORRECCIÓN 1) ---
      removeModal()
      # Usar capture.output para evitar coerción de 'environment'
      message_to_show <- if (!is.null(error$message)) error$message else paste(capture.output(print(error)), collapse = "\n")
      showModal(modalDialog(
        title = "Error en la Búsqueda NCBI",
        paste("Ocurrió un error:", message_to_show)
      ))
    })
  })
  
  
  # --- Lógica del Botón: Paso 1B (Carga de Archivo) ---
  observeEvent(input$run_file_upload, {
    
    req(input$fasta_file_upload)
    showModal(modalDialog("Cargando y validando archivo FASTA...", title = "Procesando", footer = NULL))
    
    tryCatch({
      seq_type <- input$seqType_upload
      seqs <- if (seq_type == "AA") {
        Biostrings::readAAStringSet(input$fasta_file_upload$datapath)
      } else {
        Biostrings::readDNAStringSet(input$fasta_file_upload$datapath)
      }
      
      if (length(seqs) < 3) {
        stop("El archivo FASTA debe contener al menos 3 secuencias.")
      }
      
      # 1. Guardar estado
      pipeline_state$raw_sequences <- seqs
      pipeline_state$sequence_type <- seq_type
      removeModal()
      
      # 2. Avanzar al usuario
      shinyjs::enable(selector = "#wizardTabs li a[data-value='tab_2_msa']")
      updateTabsetPanel(session, "wizardTabs", selected = "tab_2_msa")
      
    }, error = function(e) {
      # --- Error (CORRECCIÓN 2) ---
      removeModal()
      # En tryCatch, el error es 'e' y se accede con e$message
      message_to_show <- e$message
      showModal(modalDialog(
        title = "Error al Cargar Archivo",
        paste("Ocurrió un error:", message_to_show)
      ))
    })
  })
  
  
  # --- Lógica del Botón: Paso 2 (Ejecutar MSA) ---
  observeEvent(input$run_msa, {
    
    req(pipeline_state$raw_sequences)
    showModal(modalDialog(paste("Ejecutando", input$msaMethod, "..."), 
                          title = "Procesando MSA", footer = NULL))
    
    seqs <- pipeline_state$raw_sequences
    mMethod <- input$msaMethod
    
    # --- Ejecución Asíncrona ---
    future({
      AutoPhyloR::run_msa_analysis(sequences_xstringset = seqs, 
                                   msa_method = mMethod)
      # SINTAXIS CORREGIDA A (function(...) { ... })
    }) %...>% (function(result_msa) {
      # --- Éxito ---
      removeModal()
      pipeline_state$msa_alignment <- result_msa
      shinyjs::enable(selector = "#wizardTabs li a[data-value='tab_3_tree']")
      updateTabsetPanel(session, "wizardTabs", selected = "tab_3_tree")
      
    }) %...>% (function(error) {
      # --- Error (CORRECCIÓN 3) ---
      removeModal()
      # Usar capture.output para evitar coerción de 'environment'
      message_to_show <- if (!is.null(error$message)) error$message else paste(capture.output(print(error)), collapse = "\n")
      showModal(modalDialog(
        title = "Error en el MSA",
        paste("Ocurrió un error:", message_to_show)
      ))
    })
  })
  
  # Salida de texto para el MSA (Pestaña 2)
  output$msa_print_output <- renderPrint({
    req(pipeline_state$raw_sequences)
    print(pipeline_state$raw_sequences)
  })
  
  
  # --- Lógica del Botón: Paso 3 (Generar Árbol) ---
  observeEvent(input$run_tree, {
    
    req(pipeline_state$msa_alignment)
    showModal(modalDialog(paste("Generando árbol con", input$phyloMethod, "..."), 
                          title = "Procesando Inferencia", footer = NULL))
    
    msa_obj <- pipeline_state$msa_alignment
    pMethod <- input$phyloMethod
    sType <- pipeline_state$sequence_type
    
    # --- Ejecución Asíncrona ---
    future({
      AutoPhyloR::run_tree_inference(msa_alignment = msa_obj, 
                                     phylo_method = pMethod,
                                     seqs_type = sType)
      # SINTAXIS CORREGIDA A (function(...) { ... })
    }) %...>% (function(result_tree) {
      # --- Éxito ---
      removeModal()
      pipeline_state$phylo_tree <- result_tree
      shinyjs::enable(selector = "#wizardTabs li a[data-value='tab_4_results']")
      updateTabsetPanel(session, "wizardTabs", selected = "tab_4_results")
      
    }) %...>% (function(error) {
      # --- Error (CORRECCIÓN 4) ---
      removeModal()
      # Usar capture.output para evitar coerción de 'environment'
      message_to_show <- if (!is.null(error$message)) error$message else paste(capture.output(print(error)), collapse = "\n")
      showModal(modalDialog(
        title = "Error en la Inferencia",
        paste("Ocurrió un error:", message_to_show)
      ))
    })
  })
  
  
  # --- Lógica del Paso 4: Salidas y Descargas ---
  
  reactive_tree_plot_obj <- reactive({
    req(pipeline_state$phylo_tree)
    AutoPhyloR::render_tree_plot(pipeline_state$phylo_tree)
  })
  
  # 1. Renderizar el gráfico
  output$tree_plot_output <- renderPlot({
    reactive_tree_plot_obj()
  })
  
  # 2. Resumen del MSA
  output$msa_final_summary <- renderPrint({
    req(pipeline_state$msa_alignment)
    print(pipeline_state$msa_alignment)
  })
  
  # 3. Descarga de Imagen (PNG)
  output$download_plot <- downloadHandler(
    filename = function() { "AutoPhyloR_tree_plot.png" },
    content = function(file) {
      ggplot2::ggsave(file, plot = reactive_tree_plot_obj(), device = "png", 
                      width = 8, height = 10, units = "in")
    }
  )
  
  # 4. Descarga del Árbol (Newick)
  output$download_tree <- downloadHandler(
    filename = function() { "AutoPhyloR_tree.nwk" },
    content = function(file) {
      req(pipeline_state$phylo_tree)
      ape::write.tree(pipeline_state$phylo_tree, file)
    }
  )
  
  # 5. Descarga del Alineamiento (FASTA)
  output$download_msa <- downloadHandler(
    filename = function() { "AutoPhyloR_alignment.fasta" },
    content = function(file) {
      req(pipeline_state$msa_alignment)
      alignment_set <- as(pipeline_state$msa_alignment, "XStringSet")
      Biostrings::writeXStringSet(alignment_set, file)
    }
  )
}

# --- Ejecutar la Aplicación ---
shinyApp(ui = ui, server = server)