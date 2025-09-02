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
                    textInput("roi_code","Add code to subset or manipulate GRanges object"),
                    selectInput("start_feature", "Start feature", choices = c("TSS", "TES", "Exon")),
                    numericInput("start_flank", "Start flank (bp)", 0),
                    selectInput("start_direction", "Start direction", choices = c("Upstream", "Downstream")),
                    selectInput("end_feature", "End feature", choices = c("TSS", "TES", "Exon")),
                    numericInput("end_flank", "End flank (bp)", 0),
                    selectInput("end_direction", "End direction", choices = c("Upstream", "Downstream")),
                    actionButton("apply_getFeature", "Apply getFeature")
                )
              ),
              fluidRow(
                box(width = 12, title = "Regions Table", DTOutput("roi_table"),
                    downloadButton("roi_download", "Download"))
              )
      ),

      # 2. ROI Lists
      tabItem(tabName = "roi_list",
              fluidRow(
                box(width = 4, title = "Create ROI Sets",
                    uiOutput("roi_select_ui"),
                    textInput("roi_list_name", "Name for new ROI set"),
                    actionButton("save_roi_list", "Save ROI Set")
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
                box(width = 6, title = "Upload BigWig Files",
                    fileInput("bw_files", "Upload BigWigs", multiple = TRUE),
                    textInput("bw_names", "Names (comma separated)")
                ),
                box(width = 6, title = "Options",
                    checkboxInput("bw_stranded", "Stranded experiment?", value = FALSE),
                    radioButtons("strand_type", "Strand type", choices = c("rev", "for", "no"))
                )
              )
      ),

      # 4. Matrices
      tabItem(tabName = "matrix",
              fluidRow(
                box(width = 4, title = "Generate Matrices",
                    textAreaInput("extra_opts", "Additional normalizeToMatrix options", ""),
                    actionButton("gen_matrix", "Generate Matrices")
                ),
                box(width = 8, title = "Log2 Ratios",
                    selectInput("mat1", "Sample 1", choices = NULL),
                    selectInput("mat2", "Sample 2", choices = NULL),
                    numericInput("pseudo", "Pseudo count", 0.1),
                    actionButton("add_log2", "Add log2 ratio")
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

  output$roi_server <- renderUI({
    choices = list.files(input$roi_server_path,pattern = "gtf$")
    tagList(
      selectInput("roi_server_file","Select file",choices = choices)
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
    df <- eval(parse(input$roi_code))

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
      featureStart = input$start_feature,
      featureEnd = input$end_feature,
      flankStart = input$start_flank,
      flankEnd = input$end_flank
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
}

shinyApp(ui, server)
