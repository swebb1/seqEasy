# app.R
library(shiny)
library(shinydashboard)
library(DT)
library(GenomicRanges)
library(rtracklayer)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(ggplot2)
library(dplyr)
library(stringr)
library(furrr)
library(purrr)

# --- Source helper functions (you provided) ---
source("../../R/getFeature.R")
source("../../R/importBWlist.R")
source("../../R/matList.R")
source("../../R/hmList.R")
source("../../R/mplot.R")

# ------------------ UI ------------------
ui <- dashboardPage(
  dashboardHeader(title = "Genomic ROI Explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Regions of Interest", tabName = "roi", icon = icon("dna")),
      menuItem("2. ROI Lists", tabName = "roi_list", icon = icon("list")),
      menuItem("3. Import BigWigs", tabName = "bigwig", icon = icon("file-import")),
      menuItem("4. Generate Matrices", tabName = "matrix", icon = icon("th")),
      menuItem("5. Heatmaps", tabName = "heatmap", icon = icon("th-large")),
      menuItem("6. Meta Plots", tabName = "metaplot", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      # 1. ROI
      tabItem(tabName = "roi",
              fluidRow(
                box(width = 6, title = "Import Regions", status = "primary",
                    selectInput("roi_select", "From", choices = c("Server","Upload")),
                    conditionalPanel("input.roi_select == 'Server'",
                                     textInput("roi_server_path", "File path", "test_data"),
                                     uiOutput("roi_server")
                    ),
                    conditionalPanel("input.roi_select == 'Upload'",
                                     fileInput("roi_file", "Upload BED/GFF/GTF/TSV")
                    ),
                    actionButton("roi_load", "Load ROI")
                ),
                box(width = 6, title = "Filter or manipulate", status = "info",
                    textInput("roi_code","Add code to subset or manipulate GRanges object",placeholder = "df |> filter(seqnames == 'chr1')"),
                    actionButton("roi_code_apply","Apply code"),
                    selectInput("start_feature", "Start feature", choices = c("TSS", "TES", "Exon")),
                    conditionalPanel("input.start_feature == 'Exon'",
                                     numericInput("start_exon", "Start exon", 1),
                                     selectInput("start_splice","Start exon ss",choices=c("5prime","3prime"))
                    ),
                    numericInput("start_flank", "Start flank (bp)", 0),
                    selectInput("start_direction", "Start direction", choices = c("up", "down"), selected = "up"),
                    selectInput("end_feature", "End feature", choices = c("TSS", "TES", "Exon"), selected = "TES"),
                    conditionalPanel("input.end_feature == 'Exon'",
                                     numericInput("end_exon", "End exon", 1),
                                     selectInput("end_splice","End exon ss",choices=c("5prime","3prime"))
                    ),
                    numericInput("end_flank", "End flank (bp)", 0),
                    selectInput("end_direction", "End direction", choices = c("up", "down"), selected = "down"),
                    actionButton("apply_getFeature", "Apply getFeature")
                )
              ),
              fluidRow(
                box(width = 12, title = "Regions Table",
                    DTOutput("roi_table"),
                    downloadButton("roi_download", "Download"),
                    textInput("roi_list_name", "Name for new ROI set"),
                    actionButton("save_roi_list", "Save ROI Set")
                )
              )
      ),

      # 2. ROI Lists
      tabItem(tabName = "roi_list",
              fluidRow(
                box(width = 4, title = "Create ROI Sets",
                    uiOutput("roi_select_ui"),
                ),
                box(width = 8, title = "Saved ROI Sets",
                    DTOutput("roi_list_table"))
              ),
              fluidRow(
                box(width = 12, title = "Flanking Regions (single set only)",
                    conditionalPanel("output.single_set_selected == true",
                                     numericInput("flank_up", "Upstream flank (bp)", 1000),
                                     numericInput("flank_down", "Downstream flank (bp)", 1000)
                    )
                )
              )
      ),

      # 3. BigWigs
      tabItem(tabName = "bigwig",
              fluidRow(
                box(width = 6, title = "Import BigWig Files",
                    selectInput("bw_stranded", "Stranded experiment?", choices = c("Yes","No")),
                    selectInput("bw_select", "From", choices = c("Server","Upload")),
                    conditionalPanel("input.bw_select == 'Server'",
                                     textInput("bw_server_path", "File path", "test_data"),
                                     uiOutput("bw_server")
                    ),
                    conditionalPanel("input.bw_select == 'Upload'",
                                     fileInput("bw_files", "Upload BigWigs", multiple = TRUE)
                    ),
                    conditionalPanel("input.bw_stranded == 'Yes'",
                                     conditionalPanel("input.bw_select == 'Server'",
                                                      uiOutput("bw_server_rev")
                                     ),
                                     conditionalPanel("input.bw_select == 'Upload'",
                                                      fileInput("bw_files_rev", "Upload Reverse BigWigs", multiple = TRUE)
                                     ),
                    ),
                    textInput("bw_names", "Names (comma separated)"),
                    textInput("bw_names_rm","Remove from names"),
                    actionButton("load_bw", "Load BigWigs")
                ),
                box(width = 6, title = "BigWig Files",
                    DTOutput("bw_table")
                )
              )
      ),

      # 4. Matrices
      tabItem(tabName = "matrix",
              fluidRow(
                box(width = 6, title = "Generate Matrices",
                    uiOutput("matrixUI"),
                    textAreaInput("matrix_opts", "Additional normalizeToMatrix options", ""),
                    textInput("matrix_name", "Matrix Name", value = "matrix_list"),
                    actionButton("gen_matrix", "Generate Matrices")
                ),
                box(width = 6, title = "Log2 Ratios",
                    selectInput("mat1", "Sample 1", choices = NULL),
                    selectInput("mat2", "Sample 2", choices = NULL),
                    numericInput("pseudo", "Pseudo count", 0.1),
                    actionButton("add_log2", "Add log2 ratio")
                ),
                box(width = 4, title = "Matrices",
                    DTOutput("matrix_table")
                )
              )
      ),

      # 5. Heatmaps
      tabItem(tabName = "heatmap",
              fluidRow(
                box(width = 4, title = "Heatmap Options",
                    sliderInput("max_q", "Max quantile", min = 0.8, max = 1, value = 0.99),
                    sliderInput("min_q", "Min quantile", min = 0, max = 0.2, value = 0),
                    selectInput("col_fun", "Color function", choices = c("red","bl2rd","red0")),
                    checkboxInput("log2mat", "Apply log2 transform", FALSE),
                    actionButton("draw_heatmaps", "Draw Heatmaps"),
                    downloadButton("download_heatmap", "Download Heatmap")
                ),
                box(width = 8, title = "Heatmap Viewer",
                    plotOutput("heatmap_plot"),
                    InteractiveComplexHeatmapOutput("interactive_hm")
                )
              )
      ),

      # 6. Meta Plots
      tabItem(tabName = "metaplot",
              fluidRow(
                box(width = 4, title = "Metaplot Options",
                    textInput("feature_label", "Feature label", "Gene"),
                    checkboxInput("compare_control", "Compare to control?", FALSE),
                    actionButton("draw_metaplot", "Draw Metaplot"),
                    downloadButton("download_metaplot", "Download Metaplot")
                ),
                box(width = 8, title = "Metaplot",
                    plotOutput("metaplot_plot"),
                    textAreaInput("extra_layers", "Add ggplot layers", "p + theme_minimal()")
                )
              )
      )
  )
  )
)


# ------------------ SERVER ------------------
server <- function(input, output, session) {

  # --- Reactive store ---
  roi_data <- reactiveVal(NULL)   # holds imported regions (GRanges)
  roi_table_data <- reactiveVal(NULL)  # holds table version for DT
  roi_sets <- reactiveValues(list = list())   # saved ROI sets
  bwf <- reactiveValues(list = list())   # BigWig files
  bwr <- reactiveValues(list = list())   # BigWig reverse files
  matList_sets <- reactiveValues(list = list())

  output$roi_server <- renderUI({
    choices = list.files(input$roi_server_path,pattern = "gtf$")
    tagList(
      selectInput("roi_server_file","Select file",choices = choices)
    )
  })

  output$bw_server <- renderUI({
    choices = list.files(input$bw_server_path,pattern = "bw$")
    tagList(
      selectInput("bw_server_files","Select bigWig files",choices = choices,multiple = T,selectize = T)
    )
  })

  output$bw_server_rev <- renderUI({
    choices = list.files(input$bw_server_path,pattern = "bw$")
    tagList(
      selectInput("bw_server_rev_files","Select reverse bigWig files",choices = choices,multiple = T,selectize = T)
    )
  })

  output$matrixUI <- renderUI({
    tagList(
      selectInput("matrix_type","Plot type",choices = c("Feature","Meta","Combined")),
      conditionalPanel("input.matrix_type == 'Combined'",
                       selectInput("matrix_grl_combined","Select ROI list", choices = names(roi_sets$list),multiple = T,selectize = T),
      ),
      conditionalPanel("input.matrix_type != 'Combined'",
                       selectInput("matrix_grl","Select ROI list", choices = names(roi_sets$list),multiple = F,selectize = T),
                       numericInput("matrix_upstream","Extend upstream (bp)",min = 0,value = 100),
                       numericInput("matrix_downstream","Extend downstream (bp)",min = 0,value = 500),
                       numericInput("matrix_w","Window size for flanks (bp)",min = 0,value = 10),
                       selectInput("matrix_target","Include target region?", choices = c("T","F")),
                       conditionalPanel("input.matrix_target == 'T",
                                        numericInput("matrix_ratio","Target ratio",min = 0,value = 0.25)
                       )
      ),
      numericInput("matrix_wins","Number of windows per feature",min = 0,value = 10),
      selectInput("matrix_bw","Select samples", choices = names(bwf$list),multiple = T,selectize = T),
      selectInput("matrix_strand","Stranded",choices = c("no","for","rev")),
      selectInput("matrix_mode","Mode",choices = c("coverage","w0","weighted","absolute")),
      selectInput("matrix_smooth","Smooth",choices = c("T","F")),
    )
  })

  # --- Step 1: ROI import ---
  observeEvent(input$roi_load, {

    if (input$roi_select == "Server" && input$roi_server_file != "") {
      path <- file.path(input$roi_server_path,input$roi_server_file)
    } else {
      path <- input$roi_file$datapath
    }

    # Auto-detect format
    ext <- tools::file_ext(path)
    gr <- NULL
    if (ext %in% c("bed", "tsv", "txt")) {
      df <- read.delim(path, header = FALSE)
      colnames(df)[1:3] <- c("chr","start","end")
      gr <- GenomicRanges::GRanges(df)
    } else if (ext %in% c("gtf", "gff")) {
      gr <- rtracklayer::import(path)
    } else {
      showNotification("Unsupported format", type = "error")
      return(NULL)
    }

    roi_data(gr)
    roi_table_data(as.data.frame(gr))
  })

  # --- Step 1b: ROI code ---

  observeEvent(input$roi_code_apply, {
    req(input$roi_code, roi_table_data())
    df <- roi_table_data()
    df <- eval(parse(text = input$roi_code))

    gr <- GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns = T)

    roi_data(gr)
    roi_table_data(as.data.frame(gr))
  })

  # --- Step 1c: Apply getFeature ---
  observeEvent(input$apply_getFeature, {
    req(roi_data())
    gr <- roi_data()

    # Call your function
    gr_new <- getFeature(
      gr,
      start_feature = input$start_feature,
      end_feature = input$end_feature,
      start_flank = input$start_flank,
      end_flank = input$end_flank,
      start_direction = input$start_direction,
      end_direction = input$end_direction
    )

    roi_data(gr_new)
    roi_table_data(as.data.frame(gr_new))
  })

  # --- Step 1d: Show in table ---
  output$roi_table <- renderDT({
    req(roi_table_data())
    datatable(roi_table_data(), filter = "top", options = list(scrollX = TRUE))
  })

  # --- Step 1e: Download ---
  output$roi_download <- downloadHandler(
    filename = function() {
      paste0("roi_export_", Sys.Date(), ".bed")
    },
    content = function(file) {
      df <- roi_table_data()
      write.table(df, file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
  )

  # UI: list current ROI set choices
  output$roi_select_ui <- renderUI({
    sets <- names(roi_sets$list)
    selectInput("roi_set_select", "Select ROI sets", choices = sets, multiple = TRUE)
  })

  # Save current ROI into sets
  observeEvent(input$save_roi_list, {
    req(roi_data(), input$roi_list_name)
    nm <- input$roi_list_name
    roi_sets$list[[nm]] <- roi_data()
    showNotification(paste("ROI set saved:", nm), type = "message")
  })

  # Show ROI sets table
  output$roi_list_table <- renderDT({
    sets <- roi_sets$list
    if (length(sets) == 0) return(NULL)
    df <- data.frame(
      Name = names(sets),
      Count = sapply(sets, length),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t"))
  })

  # Signal if single set is selected
  output$single_set_selected <- reactive({
    length(input$roi_set_select) == 1
  })
  outputOptions(output, "single_set_selected", suspendWhenHidden = FALSE)

  # Apply flanking if single set selected
  observeEvent(c(input$flank_up, input$flank_down), {
    req(input$roi_set_select, length(input$roi_set_select) == 1)
    nm <- input$roi_set_select
    gr <- roi_sets$list[[nm]]
    gr_flanked <- GenomicRanges::resize(gr, width(gr) + input$flank_up + input$flank_down,
                                        fix = "center")
    roi_sets$list[[nm]] <- gr_flanked
    showNotification(paste("Flanking applied to set:", nm), type = "message")
  })

  # ------------------- STEP 3 (BigWig import) -------------------

  # Import BigWig file
  observeEvent(input$load_bw, {

    if (input$bw_select == "Server" && input$bw_server_path != "") {
      path <- file.path(input$bw_server_path,input$bw_server_files)
      if (input$bw_stranded == "Yes"){
        rev_path <- file.path(input$bw_server_path,input$bw_server_rev_files)
      }
    } else {
      path <- input$bw_files$datapath
      if (input$bw_stranded == "Yes"){
        rev_path <- input$bw_files_rev$datapath
      }
    }

    if (input$bw_names == ""){
      bw_names <- basename(path)
    }
    else{
      bw_names <- strsplit(input$bw_names,",")[[1]]
    }
    if (input$bw_names_rm != ""){
      bw_names <- bw_names |> str_remove(input$bw_names_rm)
    }

    bwf_import <- importBWlist(bwf = path, names = bw_names)
    bwf$list <- (bwf_import)

    if (input$bw_stranded == "Yes"){
      bwr_import <- importBWlist(bwf = rev_path, names = bw_names, selection = roi_data())
      bwr$list <- bwr_import
    }

    showNotification("BigWigs imported", type = "message")
  })

  # Show imported BigWigs
  output$bw_table <- renderDT({
    bf <- bwf$list
    br <- bwf$list
    if (length(bf) == 0) return(NULL)
    df <- data.frame(
      Forward = names(bf),
      Reverse = names(br),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t", scrollX = TRUE))
  })

  # ------------------- STEP 4 (Matrix: call matList) -------------------

  # generate matrix list
  observeEvent(input$gen_matrix, {

    grl <- roi_sets$list[[input$matrix_grl]]
    names(grl) <- input$matrix_grl

    name <- input$matrix_name

    if(input$matrix_type == "Combine"){

      grl <- roi_sets$list[[input$matrix_grl_combined]]
      names(grl) <- input$matrix_grl_combined

      wins <- str_split(input$matrix_wins,",") |> unlist() |> as.list()
      names(wins) <- names(grl)
      ml <- NULL
    }
    if(input$matrix_type == "Meta"){

      grl <- roi_sets$list[input$matrix_grl]

      extend = c(input$matrix_upstream,input$matrix_downstream)
      wins = list((sum(extend) / input$matrix_w) / (1 - input$matrix_ratio))
      names(wins) = names(grl)

      ml <- matList(bwf = bwf$list[input$matrix_bw],
                    bwr = bwr$list[input$matrix_bw],
                    names = input$matrix_bw,
                    grl = grl,
                    wins = wins,
                    mode = input$matrix_mode,
                    strand = input$matrix_strand,
                    extend = extend,
                    smooth = input$matrix_smooth |> as.logical(),
                    w = input$matrix_w,
                    include_target = input$matrix_target |> as.logical(),
                    target_ratio = input$matrix_ratio)
    }

    matList_sets$list[[name]] <- ml
    showNotification(paste("Matrix generated:",name), type = "message")

  })

  # Show generated Matrices
  output$matrix_table <- renderDT({
    mls <- matList_sets$list
    if (length(mls) == 0) return(NULL)
    df <- data.frame(
      Name = names(mls),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t", scrollX = TRUE))
  })

}

shinyApp(ui, server)
