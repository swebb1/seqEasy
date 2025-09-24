# app.R
library(shiny)
library(shinydashboard)
library(DT)
library(GenomicRanges)
library(rtracklayer)
library(genomation)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(furrr)
library(purrr)
library(dplyr)
library(readr)
library(tibble)

# --- Source helper functions (you provided) ---
source("../../R/getFeature.R")
source("../../R/importBWlist.R")
source("../../R/matList.R")
source("../../R/hmList.R")
source("../../R/mplot.R")

# ------------------ UI ------------------
ui <- dashboardPage(
  dashboardHeader(title = "SeqEasy"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Regions of Interest", tabName = "roi", icon = icon("dna")),
      menuItem("2. ROI Lists", tabName = "roi_list", icon = icon("list")),
      menuItem("3. Import BigWigs", tabName = "bigwig", icon = icon("file-import")),
      menuItem("4. Matrices", tabName = "matrix", icon = icon("th")),
      menuItem("5. Annotation", tabName = "annotation", icon = icon("pen")),
      menuItem("6. Heatmaps", tabName = "heatmap", icon = icon("th-large")),
      menuItem("7. Meta Plots", tabName = "metaplot", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      # 1. ROI
      tabItem(tabName = "roi",
              fluidRow(
                box(width = 6, title = "Import Regions", status = "primary",
                    selectInput("roi_select", "From", choices = c("Directory","Upload")),
                    conditionalPanel("input.roi_select == 'Directory'",
                                     textInput("roi_server_path", "Path", "test_data"),
                                     uiOutput("roi_server")
                    ),
                    conditionalPanel("input.roi_select == 'Upload'",
                                     fileInput("roi_file", "Upload BED/GFF/GTF/TSV")
                    ),
                    actionButton("roi_load", "Load ROI")
                ),
                box(width = 6, title = "Filter or manipulate", status = "info",
                    textInput("roi_code","Add code to subset or manipulate GRanges object",placeholder = "filter(seqnames == 'chr1')"),
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
                    textInput("roi_list_name", "Name for new ROI set"),
                    actionButton("save_roi_list", "Create ROI Set")
                )
              )
      ),

      # 2. ROI Lists
      tabItem(tabName = "roi_list",
              fluidRow(
                box(width = 4, title = "ROI Sets",
                    DTOutput("roi_list_table")
                ),
                box(width = 8, title = "ROI Set",
                    DTOutput("roi_select_table"),
                )
              ),
              fluidRow(
                box(width = 4, title = "Save ROI Set",
                    textInput("roi_select_path", "Save Path for new ROI set",value = "test_data"),
                    textInput("roi_select_name", "Save Name for new ROI set"),
                    actionButton("save_roi_select", "Save ROI Set as RDS"),
                    downloadButton("roi_download", "Download to computer")
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
                                     selectInput("bw_type", "Type", choices=c("BigWig","RDS object")),
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
                    DTOutput("bw_table"),
                    textInput("bw_save_path","Save Directory","test_data"),
                    textInput("bw_save_name", "Save Name"),
                    actionButton("save_bw", "Save BigWigs as RDS object")
                )
              )
      ),

      # 4. Matrices
      tabItem(tabName = "matrix",
              fluidRow(
                box(width = 6, title = "Generate Matrices",
                    selectInput("matrix_select","Create matrix from:", choices = c("Generate","Load")),
                    conditionalPanel("input.matrix_select == 'Load'",
                                     textInput("matrix_path", "Matrix Directory", value = "test_data"),

                    ),
                    uiOutput("matrixUI")
                ),
                box(width = 6, title = "Matrices",
                    DTOutput("matrix_table"),
                    textInput("matrix_select_path","Save Directory","test_data"),
                    textInput("matrix_select_name", "Save Name"),
                    actionButton("save_matrix_select", "Save Matrices as RDS object"),
                    div(
                      style = "max-height:300px; overflow-y:auto; background-color:#f8f9fa; padding:5px; border:1px solid #ddd;",
                      verbatimTextOutput("matrix_select_text")
                    )
                )
              ),
              fluidRow(
                box(width = 6, title = "Log2 Ratios",
                    selectInput("mat1", "Sample 1", choices = NULL),
                    selectInput("mat2", "Sample 2", choices = NULL),
                    numericInput("pseudo", "Pseudo count", 0.1),
                    actionButton("add_log2", "Add log2 ratio")
                )
              )
      ),

      # 5. Annotations
      tabItem(tabName = "annotation",
            fluidRow(
              box(width = 6, title = "Load Annotation",
                  selectInput("anno_select", "From", choices = c("Directory","Upload")),
                  conditionalPanel("input.anno_select == 'Directory'",
                                   textInput("anno_server_path", "Path", "test_data"),
                                   uiOutput("anno_serverUI"),
                                   uiOutput("anno_select_colsUI")
                  ),
                  conditionalPanel("input.anno_select == 'Upload'",
                                   fileInput("anno_file", "Upload Dataframe")
                  ),
                  textInput("anno_name", "Annotation Name",value = "Annotation"),
                  actionButton("anno_load", "Load annotation file")
              ),
              box(width = 6, title = "Annotation Table",
                  DTOutput("anno_list_table"),
                  DTOutput("anno_select_table"),
              )
            ),
            fluidRow(
              box(width = 6, title = "Annotation Colour Table",
                  DTOutput("anno_select_cols_table")
              ),
              box(width = 6, title = "Annotation Colours",
                  plotOutput("anno_select_cols_plot")
              )
            )
      ),

      # 6. Heatmaps
      tabItem(tabName = "heatmap",
              fluidRow(
                box(width = 3, title = "Heatmap Options",
                    uiOutput("heatmap_selectUI"),
                    uiOutput("heatmapUI"),
                    uiOutput("hm_annoUI"),
                    actionButton("gen_heatmaps", "Generate Heatmaps")
                ),
                box(width = 9, title = "Heatmap Viewer",
                    uiOutput("heatmap_drawUI"),
                    actionButton("draw_heatmaps","Draw Heatmaps"),
                    div(
                      style = "max-height:1000px; overflow-x:auto; background-color:#f8f9fa; padding:5px; border:1px solid #ddd;",
                      plotOutput("heatmap_plot")
                    ),
                    downloadButton("download_heatmap", "Download Heatmap")
                )
              )
      ),

      # 7. Meta Plots
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
  anno_sets <- reactiveValues(list = list())
  anno_col_sets <- reactiveValues(list = list())
  hml_data <- reactiveVal(list())
  hml_sub_data <- reactiveVal(list())

  output$roi_server <- renderUI({
    choices = list.files(input$roi_server_path,
                         pattern = "\\.(gtf|gff|tsv|txt|bed|rds)$",
                         ignore.case = TRUE)
    tagList(
      selectInput("roi_server_file","Select file",choices = choices)
    )
  })

  output$bw_server <- renderUI({
    if(input$bw_type == "BigWig"){
      choices = list.files(input$bw_server_path,pattern = "bw$",ignore.case = T)
      tagList(
        selectInput("bw_server_files","Select bigWig files",choices = choices,multiple = T,selectize = T)
      )
    }
    else{
      choices = list.files(input$bw_server_path,pattern = "rds$",ignore.case = T)
      tagList(
        selectInput("bw_server_rds","Select bigWig RDS object",choices = choices,multiple = F,selectize = T)
      )
    }
  })

  output$bw_server_rev <- renderUI({
    if(input$bw_type == "BigWig"){
      choices = list.files(input$bw_server_path,pattern = "bw$",ignore.case = T)
      tagList(
        selectInput("bw_server_rev_files","Select reverse bigWig files",choices = choices,multiple = T,selectize = T)
      )
    }
    else{
      choices = list.files(input$bw_server_path,pattern = "rds$",ignore.case = T)
      tagList(
        selectInput("bw_server_rev_rds","Select reverse bigWig RDS object",choices = choices,multiple = F,selectize = T)
      )
    }
  })

  output$matrixUI <- renderUI({
    if(input$matrix_select == "Generate"){
      tagList(
        selectInput("matrix_type","Plot type",choices = c("Meta","Combined")),
        conditionalPanel("input.matrix_type == 'Combined'",
                         selectInput("matrix_grl_combined","Select ROIs", choices = names(roi_sets$list),multiple = T,selectize = T),
                         textInput("matrix_wins_combined","Number of windows per feature (separated by a comma)")
        ),
        conditionalPanel("input.matrix_type != 'Combined'",
                         selectInput("matrix_grl","Select ROI", choices = names(roi_sets$list),multiple = F,selectize = T),
                         numericInput("matrix_upstream","Extend upstream (bp)",min = 0,value = 100),
                         numericInput("matrix_downstream","Extend downstream (bp)",min = 0,value = 500),
                         numericInput("matrix_w","Window size for flanks (bp)",min = 0,value = 10),
                         selectInput("matrix_target","Include target region?", choices = c("T","F")),
                         conditionalPanel("input.matrix_target == 'T",
                                          numericInput("matrix_ratio","Target ratio",min = 0,value = 0.25)
                         ),
                         numericInput("matrix_wins","Number of windows per feature",min = 0,value = 10),
        ),
        selectInput("matrix_bw","Select samples", choices = names(bwf$list),multiple = T,selectize = T),
        selectInput("matrix_strand","Stranded",choices = c("no","for","rev")),
        selectInput("matrix_mode","Mode",choices = c("coverage","w0","weighted","absolute")),
        selectInput("matrix_smooth","Smooth",choices = c("T","F")),
        textAreaInput("matrix_opts", "Additional normalizeToMatrix options", ""),
        textInput("matrix_name", "Matrix Name", value = "matrix_list"),
        actionButton("gen_matrix", "Generate Matrices")
      )}
    else{
      choices = list.files(input$matrix_path,pattern = "rds$",ignore.case = T)
      tagList(
        selectInput("matrix_rds", "Matrix RDS file",choices = choices,multiple = F),
        textInput("matrix_rds_name", "Name for matrix set",value = "matrix"),
        actionButton("load_matrix", "Load Matrices")
      )
    }
  })

  output$anno_serverUI <- renderUI({
    choices = list.files(input$anno_server_path,
                         pattern = "\\.(tsv|txt|csv)$",
                         ignore.case = TRUE)
    tagList(
      selectInput("anno_server_file","Select file",choices = choices)
    )
  })

  output$anno_select_colsUI <- renderUI({
    choices = list.files(input$anno_server_path,
                         pattern = "\\.(tsv)$",
                         ignore.case = TRUE)
    tagList(
      selectInput("anno_server_col_file","Select colour file",choices = c("default",choices))
    )
  })

  output$heatmap_selectUI <- renderUI({
    tagList(
      selectInput("hm_matl","Select matrix set",choices = names(matList_sets$list),selectize = T)
    )
  })

  output$heatmapUI <- renderUI({
    choices = names(matList_sets$list[[input$hm_matl]])
    choices = choices[!choices %in% "attributes"]
    hm_wins = matList_sets$list[[input$hm_matl]]$attributes$wins |> unlist() |> paste(collapse=",")
    hm_win_labels = matList_sets$list[[input$hm_matl]]$attributes$wins |> names() |> paste(collapse=",")
    hm_axis_labels = matList_sets$list[[input$hm_matl]]$attributes$labels
    anno_choices = names(anno_sets$list)

    tagList(
      selectInput("hm_select", "Select samples",choices =  choices, multiple = T,selectize = T),
      numericInput("hm_min_q","Minimum quantile",min = 0,max = 0.99,value=0),
      numericInput("hm_max_q","Maximum quantile",min = 0.01,max = 1,value=0.99),
      selectInput("hm_col_fun","Heatmap Colours", choices = c("red","red0","bl2red")),
      selectInput("hm_rownames","Show row names", choices = c(T,F)),
      textInput("hm_wins","Window lengths", value = hm_wins),
      textInput("hm_win_labels","Window labels", value = hm_win_labels),
      textInput("hm_min_y","Minimum y-axis",value="auto"),
      textInput("hm_max_y","Maximum y-axis",value="auto"),
      selectInput("hm_summarise", "Summarise by", choices = c("mean","median"),selected="mean"),
      textInput("hm_axis_labels","X-axis labels",value=hm_axis_labels),
      numericInput("hm_km","K-means clusters",value=1),
      selectInput("hm_log2","Apply log2", choices = c(T,F), selected=F),
      selectInput("hm_anno","Add annotation",choices = c("No",anno_choices), selected="F")
    )
  })

  output$hm_annoUI <- renderUI({
    choices = names(anno_sets$list[[input$hm_anno]])[-1]
    conditionalPanel("input.hm_anno != 'No'",
      selectInput("hm_anno_select","Select Annotations",choices = choices, multiple = T,selectize = T),
      selectInput("hm_anno_filter","Filter by annotation", choices = c("Yes","No")),
      selectInput("hm_anno_split","Split by annotation", choices = c("Yes","No")),
    )
  })

  output$heatmap_drawUI <- renderUI({
    tagList(
      selectInput("hm_draw_select", "Select samples",choices = names(hml_data()),multiple = T,selectize = T),
    )
  })

  # --- Step 1: ROI import ---
  observeEvent(input$roi_load, {

    if (input$roi_select == "Directory" && input$roi_server_file != "") {
      path <- file.path(input$roi_server_path,input$roi_server_file)
    } else {
      path <- input$roi_file$datapath
    }

    # Auto-detect format
    ext <- tools::file_ext(path)
    gr <- NULL
    if (ext %in% c("bed")) {
      gr <- genomation::readBed(path)
    } else if (ext %in% c("tsv", "txt")) {
      df <- read.delim(path, header = TRUE)
      gr <- GenomicRanges::GRanges(df)
    } else if (ext %in% c("gtf", "gff")) {
      gr <- rtracklayer::import(path)
    } else if (ext %in% c("rds", "RDS")) {
      gr <- readRDS(path)
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
    df <- eval(parse(text = paste("df |>",input$roi_code)))

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

  # Save current ROI into sets
  observeEvent(input$save_roi_list, {
    req(roi_data(), input$roi_list_name)
    nm <- input$roi_list_name
    roi_sets$list[[nm]] <- roi_data()
    showNotification(paste("ROI set created:", nm), type = "message")
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

  ## Show selected ROI
  output$roi_select_table <- renderDT({
    req(input$roi_list_table_rows_selected)

    # find which row was clicked
    row_idx <- input$roi_list_table_rows_selected
    item_name <- names(roi_sets$list)[row_idx]

    # extract that object
    df <- roi_sets$list[[item_name]] |> as.data.frame()

    datatable(df, options = list(scrollX = TRUE))
  })

  ## Download selected ROI
  output$roi_download <- downloadHandler(
    filename = function() {
      paste0(input$roi_select_name,".tsv")
    },
    content = function(file) {
      row_idx <- input$roi_list_table_rows_selected
      item_name <- names(roi_sets$list)[row_idx]
      df <- roi_sets$list[[item_name]] |> as.data.frame()
      write.table(df, file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
  )

  ## Save ROI as RDS
  observeEvent(input$save_roi_select,{
    req(input$roi_list_table_rows_selected, input$roi_select_name)

    row_idx <- input$roi_list_table_rows_selected
    item_name <- names(roi_sets$list)[row_idx]

    nm <- paste0(input$roi_select_name,".rds")
    path <- input$roi_select_path
    fn <- file.path(path,nm)

    saveRDS(roi_sets$list[[item_name]],fn)
    showNotification(paste("ROI set saved as RDS object:", nm), type = "message")
  })

  # ------------------- STEP 3 (BigWig import) -------------------

  # Import BigWig file
  observeEvent(input$load_bw, {
    if(input$bw_type == "BigWig"){
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
      bwr$list <- list()

      if (input$bw_stranded == "Yes"){
        bwr_import <- importBWlist(bwf = rev_path, names = bw_names, selection = roi_data())
        bwr$list <- bwr_import
      }
    }
    else{
      bwf$list <- readRDS(file.path(input$bw_server_path,input$bw_server_rds))
      bwr$list <- list()
      if(input$bw_stranded =="Yes"){
        bwr$list <- readRDS(file.path(input$bw_server_path,input$bw_server_rev_rds))
      }
    }

    showNotification("BigWigs imported", type = "message")
  })

  # Show imported BigWigs
  output$bw_table <- renderDT({
    bf <- bwf$list
    br <- bwr$list
    if (length(bf) == 0) return(NULL)
    if (length(br) == 0) {
      names_br = "NA"
    }
    else{
      names_br = names(br)
    }
    df <- data.frame(
      Forward = names(bf),
      Reverse = names_br,
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t", scrollX = TRUE))
  })

  ## Save ROI as RDS
  observeEvent(input$save_bw,{
    req(bwf$list)

    nm <- paste0(input$bw_save_name,".f.rds")
    path <- input$bw_save_path
    fn <- file.path(path,nm)
    saveRDS(bwf$list,fn)

    if(input$bw_stranded == "Yes"){
      nm <- paste0(input$bw_save_name,".r.rds")
      path <- input$bw_save_path
      fn <- file.path(path,nm)
      saveRDS(bwr$list,fn)
    }

    showNotification(paste("BigWigs saved as RDS object:", nm), type = "message")
  })

  # ------------------- STEP 4 (Matrix: call matList) -------------------

  # generate matrix list
  observeEvent(input$gen_matrix, {

    name <- input$matrix_name

    if(input$matrix_type == "Combined"){

      grl <- roi_sets$list[input$matrix_grl_combined]
      names(grl) <- input$matrix_grl_combined
      print(grl)

      wins <- str_split(input$matrix_wins_combined,",") |> unlist() |> as.list()
      wins <- wins |> map(~as.integer(.x)) #Convert wins to int
      names(wins) <- names(grl)
      print(wins)

      ml <- matList(bwf = bwf$list[input$matrix_bw],
                    bwr = bwr$list[input$matrix_bw],
                    names = input$matrix_bw,
                    grl = grl,
                    wins = wins,
                    mode = input$matrix_mode,
                    strand = input$matrix_strand,
                    smooth = input$matrix_smooth |> as.logical(),
                    attributes = T)
    }
    if(input$matrix_type == "Meta"){

      grl <- roi_sets$list[input$matrix_grl]

      extend = c(input$matrix_upstream,input$matrix_downstream)
      wins = list(( (sum(extend) / input$matrix_w) / (1 - input$matrix_ratio)) - (sum(extend) / input$matrix_w ))
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
                    target_ratio = input$matrix_ratio,
                    attributes = T)
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

  ## Save Matrix as RDS
  observeEvent(input$save_matrix_select,{
    req(input$matrix_table_rows_selected, input$matrix_select_name)

    row_idx <- input$matrix_table_rows_selected
    item_name <- names(matList_sets$list[row_idx])

    nm <- paste0(input$matrix_select_name,".rds")
    path <- input$matrix_select_path
    fn <- file.path(path,nm)

    saveRDS(matList_sets$list[[item_name]],fn)
    showNotification(paste("Matrix set saved as RDS object:", nm), type = "message")
  })

  output$matrix_select_text <- renderPrint({
    req(input$matrix_table_rows_selected)

    row_idx <- input$matrix_table_rows_selected
    matList_sets$list[row_idx]

  })

  ## Load matrix from RDS
  observeEvent(input$load_matrix,{
    req(input$matrix_path,input$matrix_rds,input$matrix_rds_name)

    name = input$matrix_rds_name
    matList_sets$list[[name]] <- readRDS(file.path(input$matrix_path,input$matrix_rds))

    showNotification(paste("Matrix loaded:",name), type = "message")
  })

  # ------------------- STEP 5 (Annotation: load annotations) -------------------

  observeEvent(input$anno_load, {

    if (input$anno_select == "Directory" && input$anno_server_file != "") {
      path <- file.path(input$anno_server_path,input$anno_server_file)
    } else {
      path <- input$anno_file$datapath
    }

    # Auto-detect format
    ext <- tools::file_ext(path)
    anno <- NULL
    if (ext %in% c("tsv")) {
      anno <- read_tsv(path, col_names = T)
    } else if (ext %in% c("txt")) {
      anno <- read_delim(path, col_names = T , delim = " ")
    } else if (ext %in% c("csv")) {
      anno <- read_csv(path, col_names = T)
    } else {
      showNotification("Unsupported format", type = "error")
      return(NULL)
    }

    #anno <- anno |> tibble::column_to_rownames(colnames(anno)[1])

    anno_sets$list[[input$anno_name]] <- anno

    cols <- NULL
    if(!input$anno_server_col_file=="default"){
      cols_df <- read_tsv(file.path(input$anno_server_path,input$anno_server_col_file),col_names = c("Column","Value","Colour"))

      cols <- cols_df |>
        mutate(Column = factor(Column,levels=Column |> unique())) |>
        group_by(Column) |>
        summarise(vec = list(setNames(Colour, Value)), .groups = "drop") |>
        deframe()
    }

    anno_col_sets$list[[input$anno_name]] <- cols

    showNotification(paste("Annotation loaded:",input$anno_name), type = "message")
  })

  # Show Anno sets table
  output$anno_list_table <- renderDT({
    sets <- anno_sets$list
    if (length(sets) == 0) return(NULL)
    df <- data.frame(
      Name = names(sets),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t"))
  })

  ## Show selected Anno
  output$anno_select_table <- renderDT({
    req(input$anno_list_table_rows_selected)

    # find which row was clicked
    row_idx <- input$anno_list_table_rows_selected
    item_name <- names(anno_sets$list)[row_idx]

    # extract that object
    df <- anno_sets$list[[item_name]] |> as.data.frame()

    datatable(df, options = list(scrollX = TRUE))
  })

  ## Show selected Anno colours table
  output$anno_select_cols_table <- renderDT({
    req(input$anno_list_table_rows_selected)

    # find which row was clicked
    row_idx <- input$anno_list_table_rows_selected
    item_name <- names(anno_col_sets$list)[row_idx]

    # extract that object
    al <- anno_col_sets$list[[item_name]]

    df <- imap_dfr(al, ~ tibble(
      Column = .y,
      Value  = names(.x),
      Colour = unname(.x)
    )) |>
      mutate(Column = factor(Column,levels = names(al)))

    datatable(df, options = list(scrollX = TRUE))
  })

  ## Show selected Anno Colours
  output$anno_select_cols_plot <- renderPlot({
    req(input$anno_list_table_rows_selected)

    # find which row was clicked
    row_idx <- input$anno_list_table_rows_selected
    item_name <- names(anno_sets$list)[row_idx]

    # extract that object
    al <- anno_col_sets$list[[item_name]]

    ## Convert to dataframe
    df_plot <- imap_dfr(al, ~ tibble(
      Column = .y,
      Value  = names(.x),
      Colour = unname(.x)
    )) |>
      mutate(Column = factor(Column,levels = names(al)))

    df_plot <- df_plot |>
      group_by(Column) |>
      mutate(y = row_number()) |>
      ungroup()

    ggplot(df_plot, aes(x = Column, y = -y, fill = Colour)) +
      geom_tile(color = "black") +
      geom_text(aes(label = Value), color = "black") +
      scale_fill_identity() +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

  })



  # ------------------- STEP 6 (Heatmap: draw heatmaps) -------------------

  observeEvent(input$gen_heatmaps, {

    matl = matList_sets$list[[input$hm_matl]][input$hm_select]

    wins = strsplit(input$hm_wins,",") |> unlist() |> as.numeric() |> as.list()
    names(wins) = strsplit(input$hm_win_labels,",") |> unlist() |> as.list()

    axis_labels = strsplit(input$hm_axis_labels,",") |> unlist()

    if(input$hm_min_y == "auto" | input$hm_max_y == "auto"){
      ylim = NULL
    }
    else{
      ylim = c(input$hm_min_y,input$hm_max_y)
    }

    if(input$hm_anno == "No"){
      hml <- hmList(matl = matl,
                    wins = wins,
                    win_labels = names(wins),
                    max_quantile = input$hm_max_q,
                    min_quantile = input$hm_min_q,
                    col_fun = input$hm_col_fun,
                    show_row_names = input$hm_rownames |> as.logical(),
                    ylim = ylim,
                    summarise_by = input$hm_summarise,
                    axis_labels = axis_labels,
                    row_km = input$hm_row_km,
                    log2 = input$hm_log2 |> as.logical()
              )
    }
    else{

      ## Filter by names in annotation file?
      if(input$hm_anno_filter == "Yes"){
       matl <- matl |> map(function(x){
            x[anno_sets$list[[input$hm_anno]] |> pull(name),]
       })
      }

      anno = data.frame(name = c(rownames(matl[[1]]))) |>
        left_join(anno_sets$list[[input$hm_anno]],by="name") |>
        column_to_rownames("name") |>
        select(all_of(input$hm_anno_select))

      anno_cols = anno_col_sets$list[[input$hm_anno]]
      anno_cols = anno_cols[c(input$hm_anno_select)]

      #if(input$hm_anno_filter == "No"){
      #  anno_cols = anno_cols |> map(function(x){x[["NA"]] ="grey";x})
      #}

      hml <- hmList(matl = matl,
                    wins = wins,
                    win_labels = names(wins),
                    max_quantile = input$hm_max_q,
                    min_quantile = input$hm_min_q,
                    col_fun = input$hm_col_fun,
                    show_row_names = input$hm_rownames |> as.logical(),
                    ylim = ylim,
                    summarise_by = input$hm_summarise,
                    axis_labels = axis_labels,
                    row_km = input$hm_row_km,
                    log2 = input$hm_log2 |> as.logical(),
                    anno = anno,
                    anno_cols = anno_cols
              )
    }
    hml_data(hml)
    showNotification(paste("Heatmap list generated"), type = "message")
  })

  output$heatmap_plot <- renderPlot({
    hml <- hml_sub_data()

    if(length(hml) == 0){
      NULL
    }

    l = length(hml)
    command = "hml[[1]]"
    if(l > 1){
        for(i in 2:l){
          command = paste0(command," + hml[[",i,"]]")
        }
    }
    eval(parse(text = paste("draw(",command,",show_heatmap_legend=T,merge_legend=T)")))

  })

  observeEvent(input$draw_heatmaps, {
    hml_sub_data(hml_data()[c(input$hm_draw_select)])
  })

}

shinyApp(ui, server)
