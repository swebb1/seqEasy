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
  dashboardHeader(title = "SeqEasy"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Regions of Interest", tabName = "roi", icon = icon("dna")),
      menuItem("2. ROI Lists", tabName = "roi_list", icon = icon("list")),
      menuItem("3. Import BigWigs", tabName = "bigwig", icon = icon("file-import")),
      menuItem("4. Matrices", tabName = "matrix", icon = icon("th")),
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

      # 5. Heatmaps
      tabItem(tabName = "heatmap",
              fluidRow(
                box(width = 4, title = "Heatmap Options",
                    uiOutput("heatmap_selectUI"),
                    uiOutput("heatmapUI"),
                    actionButton("gen_heatmaps", "Generate Heatmaps")
                ),
                box(width = 8, title = "Heatmap Viewer",
                    uiOutput("heatmap_drawUI"),
                    actionButton("draw_heatmaps","Draw Heatmaps"),
                    plotOutput("heatmap_plot"),
                    downloadButton("download_heatmap", "Download Heatmap")
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


  output$heatmap_selectUI <- renderUI({
    tagList(
      selectInput("hm_matl","Select matrix set",choices = names(matList_sets$list),selectize = T)
    )
  })

  output$heatmapUI <- renderUI({
    tagList(
      selectInput("hm_select", "Select samples",choices = names(matList_sets$list[[input$hm_matl]]), multiple = T,selectize = T),
      numericInput("hm_min_q","Minimum quantile",min = 0,max = 0.99,value=0),
      numericInput("hm_max_q","Maximum quantile",min = 0.01,max = 1,value=0.99),
      selectInput("hm_col_fun","Heatmap Colours", choices = c("red","red0","bl2red")),
      selectInput("hm_rownames","Show row names", choices = c(T,F)),
      textInput("hm_wins","Window lengths"),
      textInput("hm_win_labels","Window labels"),
      textInput("hm_max_y","Maximum y-axis",value="auto"),
      textInput("hm_min_y","Minimum y-axis",value="auto"),
      selectInput("hm_summarise", "Summarise by", choices = c("mean","median"),selected="mean"),
      textInput("hm_axis_labels","X-axis labels",value=""),
      numericInput("hm_km","K-means clusters",value=1),
      selectInput("hm_log2","Apply log2", choices = c(T,F), selected=F),
      selectInput("hm_split","Split by annotation",choices = c(T,F), selected=F),
      textInput("hm_split_cols","Annotation colours",value="")
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
    if (ext %in% c("bed", "tsv", "txt")) {
      df <- read.delim(path, header = FALSE)
      colnames(df)[1:3] <- c("chr","start","end")
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
    df <- eval(parse(text = paste("df <- ",input$roi_code)))

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

  # ------------------- STEP 5 (Heatmap: draw heatmaps) -------------------

  observeEvent(input$gen_heatmaps, {

    wins = strsplit(input$hm_wins,",") |> unlist() |> as.numeric() |> as.list()
    names(wins) = strsplit(input$hm_win_labels,",") |> unlist() |> as.list()

    axis_labels = strsplit(input$hm_axis_labels,",") |> unlist()

    if(input$hm_min_y == "auto" | input$hm_max_y == "auto"){
      ylim = NULL
    }
    else{
      ylim = c(input$hm_min_y,input$hm_max_y) |> as.numeric()
    }

    hml <- hmList(matl = matList_sets$list[[input$hm_matl]][input$hm_select],
                       wins = wins,
                       win_labels = names(wins),
                       #max_quantile = input$hm_max_q,
                       #min_quantile = input$hm_min_q,
                       #col_fun = input$hm_col_fun,
                       #show_row_names = input$hm_rownames |> as.logical(),
                       #ylim = ylim,
                       #summarise_by = input$hm_summarise,
                       #axis_labels = axis_labels,
                       #row_km = input$hm_row_km,
                       #log2 = input$hm_log2 |> as.logical()
                  )

      hml_data(hml)
      showNotification(paste("Heatmap list generated"), type = "message")
    })

    output$heatmap_plot <- renderPlot({
      hml <- hml_sub_data()

      if(input$hm_split){
        #hm <- hml$rowAnno + hml[[a]] + hml[[b]]
        NULL
      }
      else{
        draw(hml[[1]]+hml[[2]],show_heatmap_legend=T,merge_legend=T)
      }
    })

    observeEvent(input$draw_heatmaps, {
      hml_sub_data(hml_data()[c(input$hm_draw_select)])
    })

}

shinyApp(ui, server)
