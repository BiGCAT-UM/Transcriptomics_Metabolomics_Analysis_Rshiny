ui <- tagList(
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(4) {
                           float: right;
                           }
                           "))),
  fluidPage(
    
    # Allow warning/information messages
    useSweetAlert(),
    
    navbarPage(title = "IBD Analysis", id = "dataTypes",
               
               ################################################################
               
               # Transcriptomics
               
               ################################################################
               
               tabPanel("Transcriptomics Analysis", 
                        value = "transcriptomics",
                        
                        tabsetPanel(id = "tabs_trans",
                                    
                                    #***************************************************#
                                    # Data Upload
                                    #***************************************************#
                                    
                                    tabPanel("Data Upload", value = "upload_trans",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                               # Header
                                               h3(strong("Meta data")),
                                               
                                               h5("Upload the ",em("hmp2_metadata.csv"), " file"),
                                               
                                               # Upload meta data
                                               fileInput(inputId = "file1", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv",
                                                                    ".tsv")),
                                               
                                               # Input: Checkbox if file has header
                                               checkboxInput("header1", "Header", TRUE),
                                               
                                               # Input: Checkbox if file has row names
                                               checkboxInput("rowNames1", "Rownames", TRUE),
                                               
                                               # Input: Select separator
                                               radioButtons("sep1", "Separator",
                                                            choices = c(Comma = ",",
                                                                        Semicolon = ";",
                                                                        Tab = "\t"),
                                                            selected = ","),
                                               # Horizontal line
                                               tags$hr(),
                                               
                                               # Header
                                               h3(strong("Host transcriptomics")),
                                               
                                               h5("Upload the ",em("host_tx_counts.tsv"), " file"),
                                               
                                               # Upload file
                                               fileInput(inputId = "file2", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv",
                                                                    ".tsv")),
                                               
                                               # Input: Checkbox if file has header
                                               checkboxInput("header2", "Header", TRUE),
                                               
                                               # Input: Checkbox if file has row names
                                               checkboxInput("rowNames2", "RowNames", TRUE),
                                               
                                               # Input: Select separator
                                               radioButtons("sep2", "Separator",
                                                            choices = c(Comma = ",",
                                                                        Semicolon = ";",
                                                                        Tab = "\t"),
                                                            selected = "\t"),
                                               
                                               # Horizontal line
                                               tags$hr(),
                                               
                                               #Go forward
                                               actionBttn(inputId ="upload_NEXT", 
                                                          label ="Next", 
                                                          style = "jelly",
                                                          btn_type = "button", 
                                                          color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               # Output: Data file-1 ----
                                               uiOutput("metaText"),
                                               DT::dataTableOutput("fileContent1")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                              
                                               uiOutput("countText"),
                                               # Output: Data file-2 ----
                                               DT::dataTableOutput("fileContent2")%>% 
                                                 withSpinner(color="#0dc5c1")
                                               
                                             )
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Filtering
                                    #***************************************************#
                                    
                                    tabPanel("Filtering", value = "filtering_trans",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               # Title + description
                                               h3(strong("Filtering")),
                                               hr(),
                                               h4(strong("Sample and gene filtering")),
                                               h5("1. Samples with all zero values across all genes will be filtered."),
                                               h5("2. Genes with all zero values across all samples will be filtered."),
                                               br(),
                                               
                                               # Title + description
                                               h4(strong("Average log CPM (Count per million) filtering")),
                                               h5("This filtering step will be performed to remove lowly expressed genes."),
                                               
                                               # log CPM threshold
                                               numericInput(inputId = "threshold", 
                                                            label = "Log CPM threshold", 
                                                            value = 1),
                                               
                                               hr(),
                                               actionBttn(inputId ="cpm_filtering", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               #Go forward
                                               hr(),
                                               actionBttn(inputId ="preprocess_NEXT", 
                                                          label ="Next", 
                                                          style = "jelly",
                                                          btn_type = "button", 
                                                          color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               # Histogram
                                               uiOutput("histText_pre"),
                                               plotOutput("plot_prefiltering")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                               uiOutput("histText_post"),
                                               plotOutput("plot_postfiltering")%>% 
                                                 withSpinner(color="#0dc5c1")
                                             )
                                    ),
                                    
                                    
                                    #***************************************************#
                                    # Normalization & QC
                                    #***************************************************#
                                    
                                    tabPanel("Normalization & QC" , value = "norm_trans",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                               # Title + description
                                               h3(strong("Normalization & QC")),
                                               hr(),
                                               
                                               # Remove outliers
                                               h4(strong("Exclude Samples")),
                                               h5("Here you can exclude samples from the subsequent analysis."),
                                               awesomeCheckbox(inputId = "outlierCheckBox",
                                                               label = "Keep all samples",
                                                               value = TRUE,
                                                               status = "primary"),
                                               
                                               uiOutput("outliersout"),
                                               hr(),
                                               
                                               # Normalization and QC
                                               h4(strong("Normalization")),
                                               h5("Normalization using the DESeq package will be performed."),
                                               actionBttn(inputId ="outlierButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               hr(),
                                               
                                               #Go forward
                                               actionBttn(inputId ="norm_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput(inputId = "whichQCplot",
                                                                    label = "QC plot",
                                                                    choices = c("PCA",
                                                                                "Boxplot",
                                                                                "Heatmap"),
                                                                    selected = "PCA")
                                                 ),
                                                 column(width = 3,
                                                        radioButtons(inputId = "normalizedQC", 
                                                                     label = "Raw or Normalized Data",
                                                                     choices = c("Raw",
                                                                                 "Normalized"),
                                                                     selected = "Normalized")
                                                        ),
                                                 column(width = 3,
                                                        radioButtons(inputId = "colorQC", 
                                                                     label = "Color by",
                                                                     choices = c("Location",
                                                                                 "Disease"),
                                                                     selected = "Location")
                                                        )
                                               ),
                                               hr(),
                                               
                                               imageOutput("QCplot",
                                                           width = "700px",
                                                           height = "auto")
                                               
                                             )#End of mainPanel
                                            
                                    ),
                                    
                                    #***************************************************#
                                    # DEG analysis
                                    #***************************************************#
                                    tabPanel("DEG Analysis" , value = "deg_trans",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                               # Title + description
                                               h3(strong("Differential Gene Expression Analysis")),
                                               h5("Statistical Analysis to identify differentially expressed genes (DEGs)."),
                                               hr(),
                                               radioButtons(inputId = "AdjOrRaw", 
                                                            label = "Adjusted or raw p-value",
                                                            choices = c("Raw",
                                                                        "Adjusted"),
                                                            selected = "Raw"),
                                               br(),
                                               numericInput(
                                                 inputId = "pthreshold",
                                                 label = "P-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               br(),
                                               
                                               numericInput(
                                                 inputId = "FCthreshold",
                                                 label = "FC threshold",
                                                 value = 1.5,
                                                 min = 1,
                                                 max = 10,
                                                 step = 0.1
                                               ),
                                               br(),
                                               h5("WARNING: The selected thresholds will also be used for the 
                                                  pathway and network analysis. So, choose your thresholds carefully!"),
                                               
                                               actionBttn(inputId ="DEGButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                              
                                               hr(),
                                               actionBttn(inputId ="deg_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               
                                               selectInput(inputId = "Comparison",
                                                           label = NULL,
                                                           choices = c("Ileum: CD vs non-IBD",
                                                                       "Rectum: CD vs non-IBD",
                                                                       "Ileum: UC vs non-IBD",
                                                                       "Rectum: UC vs non-IBD"),
                                                           selected = "Ileum: CD vs non-IBD"),
                                               
                                               #output from DE analysis
                                               DT::dataTableOutput("topTable")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                               
                                               plotOutput("VolcanoPlot", width = "800px", height = "700px")%>% 
                                                 withSpinner(color="#0dc5c1")
                                                 
                                             )#End of mainPanel
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Identifier Mapping
                                    #***************************************************#
                                    
                                    tabPanel("Identifier Mapping", value = "mapping_trans",
                                             br(),
                                             #==========================================#
                                             # side panel
                                             #==========================================#
                                             sidebarPanel(
                                               #Go forward
                                               #******************************************************#
                                               #   Information
                                               #******************************************************#
                                               
                                               h3(strong("Identifier Mapping")),
                                               h5("HGNC gene symbols will be mapped to ENTREZ and ENSEMBL identifiers."),
                                               br(),
                                               
                                               actionBttn(inputId ="mappingButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),

                                               # Horizontal line ----
                                               tags$hr(),
                                               actionBttn(inputId ="mapping_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               selectInput(inputId = "mappingDisease",
                                                           label = NULL,
                                                           choices = c("CD", "UC"),
                                                           selected = "CD"),
                                               
                                               #output from DE analysis
                                               DT::dataTableOutput("mappingTable")
                                             )
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Pathway Analysis
                                    #***************************************************#
                                    
                                    tabPanel("Pathway Analysis", value = "pathway_trans",
                                             br(),
                                             sidebarPanel(
                                               #==========================================#
                                               # side panel
                                               #==========================================#
                                               h3(strong("Pathway Analysis")),
                                               h5("WikiPathways Overrepresentation Analysis will be performed. 
                                                  The logFC and p-value thresholds from the DEG Analysis will be used."),
                                               br(),
                                               
                                               actionBttn(inputId ="pathwayButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               # Horizontal line ----
                                               tags$hr(),
                                               actionBttn(inputId ="pathway_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               selectInput(inputId = "pathwayComparison",
                                                           label = NULL,
                                                           choices = c("Ileum: CD vs non-IBD",
                                                                       "Rectum: CD vs non-IBD",
                                                                       "Ileum: UC vs non-IBD",
                                                                       "Rectum: UC vs non-IBD"),
                                                           selected = "Ileum: CD vs non-IBD"),
                                               
                                               #output from DE analysis
                                               plotOutput("pathwayPlot",
                                                          width = "1000px",
                                                          height = "400px"),
                                               hr(),
                                               DT::dataTableOutput("pathwayTable")
                                             )
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Heatmap
                                    #***************************************************#
                                    
                                    tabPanel("Create Heatmap", value = "heatmap_trans",
                                             br(),
                                             sidebarPanel(
                                               #==========================================#
                                               # side panel
                                               #==========================================#
                                               h3(strong("Heatmap")),
                                               h5("Heatmap for the pathway analysis will be created"),
                                               br(),
                                               numericInput(
                                                 inputId = "p_threshold_pathway",
                                                 label = "adj.p-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               br(),
                                               numericInput(
                                                 inputId = "nSignGenesPathways",
                                                 label = "# Significant genes",
                                                 value = 5,
                                                 min = 0,
                                                 max = 20,
                                                 step = 1
                                               ),
                                               
                                               actionBttn(inputId ="heatmapButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               # Horizontal line ----
                                               tags$hr(),
                                               actionBttn(inputId ="heatmap_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               imageOutput("HeatmapPlot",
                                                           width = "700px",
                                                           height = "auto")
                                             )
                                    ),
                                    
                                    #***************************************************#
                                    # Network analysis
                                    #***************************************************#
                                    
                                    tabPanel("Network Analysis", value = "network_trans",
                                             br(),
                                             sidebarPanel(
                                               #==========================================#
                                               # side panel
                                               #==========================================#
                                               h3(strong("Network Analysis")),
                                               h5("Network analysis will be performed.  
                                                  Please make sure to have Cytoscape opened before running 
                                                  the analysis."),
                                               br(),
                                               numericInput(
                                                 inputId = "PPI_cutoff",
                                                 label = "STRING protein query cutoff",
                                                 value = 0.7,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.1
                                               ),
                                               br(),
                                               actionBttn(inputId ="networkButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary")
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               selectInput(inputId = "location_network",
                                                           label = NULL,
                                                           choices = c(Ileum = "ileum", 
                                                                       Rectum = "rectum"),
                                                           selected = "ileum"),
                                               plotlyOutput("NetworkPlot"),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               imageOutput("legendNetwork")
                                             )
                                      )
                                    
                        ) # End of tabset panel
                        
               ),#tabPanel-trans
               
               ################################################################
               
               # Metabolomics
               
               ################################################################
               tabPanel("Metabolomics Analysis", 
                        value = "metabolomics",
                        
                        tabsetPanel(id="tabs_mets",
                                    
                                    #***************************************************#
                                    # Data Upload
                                    #***************************************************#
                                    
                                    tabPanel("Data Upload", value = "upload_mets",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                               # Header
                                               h3(strong("Meta data")),
                                               
                                               h5("Upload the ",em("hmp2_metadata.csv"), " file"),
                                               
                                               # Upload meta data
                                               fileInput(inputId = "metaFile", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv")),
                                               
                                               # Input: Checkbox if file has header
                                               checkboxInput("headerMet", "Header", TRUE),
                                               
                                               # Input: Checkbox if file has row names
                                               checkboxInput("rowNamesMet", "Rownames", TRUE),
                                               
                                               # Input: Select separator
                                               radioButtons("sepMet", "Separator",
                                                            choices = c(Comma = ",",
                                                                        Semicolon = ";",
                                                                        Tab = "\t"),
                                                            selected = ","),
                                               # Horizontal line
                                               tags$hr(),
                                               
                                               
                                               
                                               h3(strong("Metabolomics data")),
                                               
                                               h5("Upload the ",em("metabolomics.csv"), " file"),
                                               
                                               # Upload file
                                               fileInput(inputId = "metsData", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv",
                                                                    ".tsv")),
                                               
                                               # Input: Checkbox if file has header
                                               checkboxInput("headMet", "Header", TRUE),
                                               
                                               # Input: Checkbox if file has row names
                                               checkboxInput("rowsMet", "RowNames", TRUE),
                                               
                                               # Input: Select separator
                                               radioButtons("sepMet2", "Separator",
                                                            choices = c(Comma = ",",
                                                                        Semicolon = ";",
                                                                        Tab = "\t"),
                                                            selected = ","),
                                               
                                               # Horizontal line
                                               tags$hr(),
                                               
                                               #Go forward
                                               actionBttn(inputId ="metUpload_NEXT", 
                                                          label ="Next", 
                                                          style = "jelly",
                                                          btn_type = "button", 
                                                          color = "danger",
                                                          icon = icon("arrow-right")
                                               )
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               # Output: Data file-1 ----
                                               uiOutput("mbxMetaText"),
                                               DT::dataTableOutput("mbxMetaFile")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                               uiOutput("mbxDataFileText"),
                                               DT::dataTableOutput("mbxDataFile")%>% 
                                                 withSpinner(color="#0dc5c1")
                                               
                                               
                                             )#mainPanel
                                             
                                    ),#tabPanel
                                    
                                    #***************************************************#
                                    # Filtering
                                    #***************************************************#
                                    tabPanel("Filtering",value = "filtering_mets",
                                             
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               # Title + description
                                               h3(strong("Filtering")),
                                               hr(),
                                               h4(strong("Sample and gene filtering")),
                                               h5("1. Samples will be filtered based on visit number and data type."),
                                               h5("2. Metabolites with all zero values across all samples will be filtered."),
                                               h5("3. Metabolites having >50% missing values will be filtered."),
                                               br(),
                                               
                                               # Apply filtering
                                               actionBttn(inputId ="metsFiltering", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               tags$hr(),
                                               
                                               
                                               #Go forward
                                               actionBttn(inputId ="metPre_NEXT", 
                                                          label ="Next", 
                                                          style = "jelly",
                                                          btn_type = "button", 
                                                          color = "danger",
                                                          icon = icon("arrow-right")
                                               )
                                               
                                             ),#sidebarPanel
                                             mainPanel(
                                               # Output: Data file-1 ----
                                               uiOutput("CDpreprocessText"),
                                               DT::dataTableOutput("mbxCDPreprocessed")%>%
                                                 withSpinner(color="#0dc5c1"),
                                               br(),
                                               # Output: Data file-2 ----
                                               uiOutput("UCpreprocessText"),
                                               DT::dataTableOutput("mbxUCPreprocessed")%>%
                                                 withSpinner(color="#0dc5c1")
                                               
                                             )#mainPanel
                                             
                                    ),#tabPanel
                                    
                                    #***************************************************#
                                    # Normalization & QC for metabolomics
                                    #***************************************************#
                                    
                                    tabPanel("Normalization" , value = "norm_mets",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                               # Title + description
                                               h3(strong("Normalization & QC")),
                                               hr(),
                                               h3(strong("Normalization" )),
                                               h5("Please select a normalization method."),
                                               selectInput(inputId = "whichNormMethod",
                                                           label = NULL,
                                                           choices = c("log2 transformation",
                                                                       "log10 transformation",
                                                                       "cube root transformation",
                                                                       "square root transformation"),
                                                           selected = "log2 transformation"),
                                               
                                               actionBttn(inputId ="normButton", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               hr(),
                                               
                                               #Go forward
                                               actionBttn(inputId ="metNorm_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               
                                               uiOutput("CDhistogram"),
                                               verbatimTextOutput("CD_shapiro_result",placeholder = TRUE),
                                               selectInput(inputId = "whichHistCD",
                                                           label = NULL,
                                                           choices = c("Normalized",
                                                                       "Raw"),
                                                           selected = "Normalized"),
                                               
                                               imageOutput("histPlotCD",
                                                           width = "800px",
                                                           height = "auto"),
                                               
                                               uiOutput("UChistogram"),
                                               verbatimTextOutput("UC_shapiro_result",placeholder = TRUE),
                                               selectInput(inputId = "whichHistUC",
                                                           label = NULL,
                                                           choices = c("Normalized",
                                                                       "Raw"),
                                                           selected = "Normalized"),
                                               
                                               imageOutput("histPlotUC",
                                                           width = "800px",
                                                           height = "auto")
                                               
                                               
                                               
                                             )#End of mainPanel
                                             
                                    ),#tabPanel-norm
                                    
                                    #***************************************************#
                                    # Statistical Analysis
                                    #***************************************************#
                                    tabPanel("Statistical Analysis", value= "stat_mets",
                                             br(),

                                      sidebarPanel(

                                        # Title + description
                                        h3(strong("Statistical Analysis")),
                                        h5("Statistical Analysis to identify statistically changed metabolites"),
                                        hr(),

                                        numericInput(
                                          inputId = "pthresholdMet",
                                          label = "P-value threshold",
                                          value = 0.05,
                                          min = 0,
                                          max = 1,
                                          step = 0.01
                                        ),
                                        br(),

                                        numericInput(
                                          inputId = "FCthresholdMet",
                                          label = "FC threshold",
                                          value = 2,
                                          min = 1,
                                          max = 10,
                                          step = 0.1
                                        ),
                                        br(),

                                        actionBttn(inputId ="statButton", label ="Apply", style = "jelly",
                                                   btn_type = "button", type = "primary", color = "primary"),

                                        hr(),
                                        actionBttn(inputId ="stat_NEXT", label ="Next", style = "jelly",
                                                   btn_type = "button", type = "primary", color = "danger",
                                                   icon = icon("arrow-right"))


                                      ),#sidebarPanel
                                      
                                      mainPanel(

                                        selectInput(inputId = "metCompPair",
                                                    label = NULL,
                                                    choices = c("CD vs non-IBD",
                                                                "UC vs non-IBD"),
                                                    selected = "CD vs non-IBD"
                                                    ),

                                        DT::dataTableOutput("metResTable")%>%
                                          withSpinner(color="#0dc5c1"),

                                        plotOutput("metVolcanoPlot", width = "800px", height = "700px")%>%
                                          withSpinner(color="#0dc5c1")

                                      )#mainPanel
                                      

                                    ),#tabPanel-stat
                                    
                                    # #***************************************************#
                                    # # Identifier mapping
                                    # #***************************************************#
                                    tabPanel("Identifier Mapping", value = "mapping_mets",
                                      
                                             br(),
                                             #==========================================#
                                             # side panel
                                             #==========================================#
                                             sidebarPanel(
                                               #Go forward
                                               #******************************************************#
                                               #   Information
                                               #******************************************************#
                                               
                                               h3(strong("Identifier Mapping")),
                                               h5("To visualize the multi-omics data in CytoScape, metabolomic data need to be prepared 
                                               for importing according to the unified database identifiers available.
                                               Therefore, metabolite data is mapped to ChEBI IDs"),
                                               br(),
                                               
                                               h4(strong("Data download")),
                                               h5("To start with mapping you should first download metabolites.bridge file  ", 
                                                  tags$a(href = "https://figshare.com/ndownloader/files/26001794", "here")  
                                                  ," and locate it to /data folder. NOTE: If you already downloaded the .bridge file, 
                                                  you can continue by clicking the Apply button! "),
                                               br(),
                                               
                                               actionBttn(inputId ="mappingButtonMets", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               # Horizontal line ----
                                               tags$hr(),
                                               actionBttn(inputId ="mappingMets_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               selectInput(inputId = "metCompPairMapping",
                                                           label = NULL,
                                                           choices = c("CD vs non-IBD",
                                                                       "UC vs non-IBD"),
                                                           selected = "CD vs non-IBD"
                                               ),
                                               
                                               #output from identifier mapping
                                               DT::dataTableOutput("metMappingTable")
                                             )         
                                             
                                    ),
                                    
                                    
                                    
                                    # #***************************************************#
                                    # # Pathway Analysis
                                    # #***************************************************#
                                      tabPanel("Pathway Analysis", value = "pathway_mets",
                                        br(),       
                                        sidebarPanel(
                                          
                                          #==========================================#
                                          # side panel
                                          #==========================================#
                                          h3(strong("Pathway Analysis")),
                                          h5("WikiPathways Overrepresentation Analysis will be performed. 
                                                  The logFC and p-value thresholds from the statistical analysis will be used."),
                                          br(),
                                          
                                          actionBttn(inputId ="pathwayButtonMet", label ="Apply", style = "jelly",
                                                     btn_type = "button", type = "primary", color = "primary"),
                                          
                                        ),#sidePanel
                                        
                                        mainPanel(
                                          selectInput(inputId = "metCompPairPathway",
                                                      label = NULL,
                                                      choices = c("CD vs non-IBD",
                                                                  "UC vs non-IBD"),
                                                      selected = "CD vs non-IBD"
                                          ),
                                          
                                          DT::dataTableOutput("metPathwayRes")%>%
                                            withSpinner(color="#0dc5c1"),
                                        )
                                      ),
                                    # 
                                 
                                    
                                    
                        )#tabSetPanel-mets
               ),#tabPanel-mets
               
               ################################################################
               #
               # # Multi-omics
               # 
               # ################################################################
               tabPanel("Multi-omics Visualization",
                        value = "multiomics",
                        tabsetPanel(id="tabs_multi",

                                    #***************************************************#
                                    # Pathway Selection
                                    #***************************************************#
                                    tabPanel("Pathway Selection",value = "selection",
                                             br(),
                                             #==========================================#
                                             # Side bar panel
                                             #==========================================#
                                             sidebarPanel(
                                               h3(strong("Pathway Selection")),
                                               h5("Select pathways for multi-omics visualization."),
                                               br(),
                                               h4(strong("Transcriptomics Thresholds")),
                                               numericInput(
                                                 inputId = "p_threshold_multi_trans",
                                                 label = "Adj. p-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               numericInput(
                                                 inputId = "q_threshold_multi_trans",
                                                 label = "q-value threshold",
                                                 value = 0.02,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               numericInput(
                                                 inputId = "nProteinsPathwayTrans",
                                                 label = "# Proteins",
                                                 value = 5,
                                                 min = 0,
                                                 max = 20,
                                                 step = 1
                                               ),
                                               hr(),
                                               h4(strong("Metabolomics Thresholds")),
                                               numericInput(
                                                 inputId = "p_threshold_multi_mets",
                                                 label = "p-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               numericInput(
                                                 inputId = "nMetsPathway",
                                                 label = "# Significant Metabolites",
                                                 value = 3,
                                                 min = 0,
                                                 max = 20,
                                                 step = 1
                                               ),
                                               numericInput(
                                                 inputId = "nProteinsPathwayMets",
                                                 label = "# Significant genes",
                                                 value = 5,
                                                 min = 0,
                                                 max = 20,
                                                 step = 1
                                               ),
                                               br(),
                                               actionBttn(inputId ="selectionButtonMulti", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               tags$hr(),
                                               actionBttn(inputId ="selectionMulti_NEXT", label ="Next", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "danger",
                                                          icon = icon("arrow-right"))
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput(inputId = "whichDisease_multi",
                                                                    label = "Disease",
                                                                    choices = c("CD",
                                                                                "UC"),
                                                                    selected = "CD")
                                                 ),
                                                 column(width = 3,
                                                        selectInput(inputId = "whichLocation_multi", 
                                                                     label = "Location",
                                                                     choices = c("Ileum",
                                                                                 "Rectum"),
                                                                     selected = "Ileum")
                                                 ),
                                               ),
                                               hr(),
                                               uiOutput("selectedPathways_header"),
                                               DT::dataTableOutput("selectedPathways_multi")%>%
                                                 withSpinner(color="#0dc5c1")
                                             )
                                             ),

                                    #***************************************************#
                                    # Visualization
                                    #***************************************************#
                                    tabPanel("Visualization", value = "visualization",
                                             br(),
                                             #==========================================#
                                             # Side bar panel
                                             #==========================================#
                                             sidebarPanel(
                                               h3(strong("Multi-omics Visualization")),
                                               br(),
                                               textInput(inputId = "vis_pathway",
                                                         label = "Pathway ID",
                                                         value = "WP15"),
                                               radioButtons(inputId = "vis_location", 
                                                            label = "Location",
                                                            choices = c(Ileum = "ileum", 
                                                                        Rectum = "rectum"),
                                                            selected = "ileum"),
                                               radioButtons(inputId = "vis_disease", 
                                                            label = "Disease",
                                                            choices = c(CD = "CD", 
                                                                        UC = "UC"),
                                                            selected = "CD"),
                                               br(),
                                               actionBttn(inputId ="visButtonMulti", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary")
                                             ),
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               plotlyOutput("multiPlot"),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               br(),
                                               plotOutput("legendMulti")
                                             )
   
                                             )

                        )#tabsetPanel
               ),#tabPanel
               
               ################################################################
               #
               # # Documentation
               # 
               # ################################################################
               tabPanel("Documentation",
                        value = "Documentation",
                        icon = icon("far fa-question-circle"),
                        tabsetPanel(id="tabs_documentation",
                                    
                                    #***************************************************#
                                    # Documentation of transcriptomics analysis
                                    #***************************************************#
                                    tabPanel("Transcriptomics Analysis",value = "doc_trans",
                                             h2(strong("Documentation of Transcriptomics Analysis")),
                                             hr(),
                                             
                                             
                                             h3(strong("1. Data pre-processing")),
                                             br(),
                                             h4(strong("1.1. Data upload")),
                                             h5("Before starting with the transcriptomics analysis,
                                                the data needs to be uploaded first. Particularly, the 
                                                transcriptomics analysis requires two files:"),
                                             h5(strong("1. Meta data file: "), em("hmp2_metadata.csv")),
                                             h5(strong("2. RNA-seq count data file: "), 
                                                em("host_tx_counts.tsv")),
                                             h5("Both files can be found in the ", em("'data'"), " folder in
                                                the working directory. The default settings for the data upload are
                                                correct and, thus, only the file location needs to be selected manually."),
                                             br(),
                                             
                                             
                                             h4(strong("1.2. Filtering")),
                                             h5("In this part, three filtering steps will be applied:"), 
                                             h5("1. Removal of samples with zero values only."),
                                             h5("2. Removal of genes with zero values only."),
                                             h5("3. Removal of genes with a mean logCPM below the selected threshold. 
                                                For this, the default threshold is 1 logCPM, but can be changed by the user."),
                                             h5("NOTE: The filtered meta data and RNA-seq count data can be found in the ",
                                                em("'1-data_preprocessing'"), " folder."),
                                             br(),
                                             
                                             
                                             h4(strong("1.3. Normalization & QC")),
                                             h5("Using the filtered data, normalization with the ", 
                                                tags$a(href = "https://bioconductor.org/packages/release/bioc/html/DESeq2.html", 
                                                        "DESeq2")  
                                                ," package will be performed. Additional quality control 
                                                will be performed through ",
                                                tags$a(href = "https://doi.org/10.1093/nar/gkt293", 
                                                       "ArrayAnalysis"), "
                                                functions (v.2). The QC plots of the raw and normalized data are saved in the ",
                                                em("'2-differential_gene_expression_analysis/QCraw'"), " and ",
                                                em("'2-differential_gene_expression_analysis/QCnorm'"), " folders, respectively."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("2. Differential gene expression analysis")),
                                             h5("Also the differential gene expression (DEG) analysis will be performed using the ", 
                                                tags$a(href = "https://bioconductor.org/packages/release/bioc/html/DESeq2.html", 
                                                       "DESeq2")  
                                                ," package. All output plots and tables can be found in the ",
                                                em("'2-differential_gene_expression_analysis/statsmodel'"), " folder."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("3. Identifier mapping")),
                                             h5("In the identifier mapping procedure, gene symbols will 
                                             be converted to Entrez Gene and ENSEMBL IDs. For this, the ",
                                             tags$a(href = "https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html", 
                                                       "AnnotationDbi"),
                                             "package (v1.58.0) will be used. In the applied methodology, 
                                             one-to-many mappings will be omitted by selecting the first hit only. 
                                             The identifier mapping output is located in the ",
                                             em("'3-identifier_mapping'"), " folder."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("4. Pathway analysis")),
                                             h5("In this step, pathway overrepresentation analysis (ORA) will be performed 
                                             using the mapped identifiers (Identifier Mapping step) and the selected logFC and 
                                             p-value thresholds (DEG Analysis step). The ORA implementation in ",
                                             tags$a(href = "https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html", 
                                                      "clusterProfiler"),
                                             " (v4.2.2) package will be applied to perform ORA 
                                             on WikiPathways, version 20220510 (747 pathways). For this analysis, 
                                             the WikiPathways .gmt file needs to be located in the ", em("'4-pathway_analysis'"), 
                                             " folder. The false discovery rate (FDR) approach will be used to acquire adjusted p-values. 
                                             Accordingly, the q-values were calculated to prevent a high FDR in multiple testing. 
                                             The output files can be found in the ", em("'4-pathway_analysis'"), " folder."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("5. Heatmap")),
                                             h5("A heatmap of pathway significance for the different 
                                                biopsy locations and diseases will be created using the ",
                                                tags$a(href = "https://github.com/raivokolde/pheatmap", 
                                                       "pheatmap"), " package."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("6. Network analysis")),
                                             h5("This part requires Cytoscape to be installed and opened. 
                                             The first step of the network analysis includes the construction 
                                             of Protein-Protein-Interaction (PPI) networks from the STRING database 
                                             using the stringApp (v2.0.0) in Cytoscape through the R package ",
                                             tags$a(href = "https://bioconductor.org/packages/release/bioc/html/RCy3.html", 
                                                       "RCy3"), " (v2.14.1). 
                                             The STRING database will be queried for the overlapping DEGs per biopsy location 
                                             with a confidence score selected by the user (default = 0.7) to construct a PPI network, 
                                             extended with the curated human pathway collection from WikiPathways, version 20220210, 2018) 
                                             to create Protein-Protein-Pathway interaction (PPP-I) networks. 
                                             The DEG data, including log2FC for CD and UC, 
                                             will be added to the networks for each biopsy location individually, 
                                             by linking Entrez gene IDs. The gene expression changes of CD and UC will be visualized 
                                             to allow for the investigation of up-regulated and down-regulated genes in both diseases. 
                                             The Markov Clustering Algorithm (MCL) was applied to the PPP-I networks through the 
                                             Cytoscape plugin named clusterMaker2 (v2.0) to investigate similarities 
                                             between genes and pathways."),
                                             br()
                                             ),
                                    
                                    #***************************************************#
                                    # Documentation of Metabolomics analysis
                                    #***************************************************#
                                    tabPanel("Metabolomics Analysis",value = "doc_mets",
                                             h2(strong("Documentation of Metabolomics Analysis")),
                                             hr(),
                                             
                                             h3(strong("7. Data pre-processing")),
                                             br(),
                                             h4(strong("7.1. Data upload")),
                                             h5("Before starting with the metabolomics analysis,
                                                the data needs to be uploaded first. Particularly, the 
                                                metabolomics analysis requires two files:"),
                                             h5(strong("1. Meta data file: "), em("hmp2_metadata.csv")),
                                             h5(strong("2. Metabolomics data file: "), 
                                                em("metabolomics.csv")),
                                             h5("Both files can be found in the ", em("'data'"), " folder in
                                                the working directory. The default settings for the data upload are
                                                correct and, thus, only the file location needs to be selected manually."),
                                             br(),
                                             
                                             
                                             h4(strong("7.2. Filtering")),
                                             h5("In this part, three filtering steps will be applied:"), 
                                             h5("1. Samples will be filtered based on visit number and data type."),
                                             h5("2. Metabolites with all zero values across all samples will be filtered."),
                                             h5("3. Metabolites having >50% missing values will be filtered."),
                                             h5("NOTE: The filtered metabolomics data can be found in the ",
                                                em("'7-metabolite_data_preprocessing/filtered'"), " folder."),
                                             br(),
                                             
                                             
                                             h4(strong("7.3. Normalization & QC")),
                                             h5("In this step, different normalization methods including cube-root, square-root, 
                                                log2, and log10 transformation can be applied. Shapiro-Wilk normality test will be 
                                                applied to observe whether the normalization method generates data with normal distribution. 
                                                The QC plots of the raw and normalized data are saved in the ",
                                                em("'7-metabolite_data_preprocessing/normalized'"), " folder."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("8. Significantly changed metabolites analysis")),
                                             h5("In the statstical analysis procedure, 
                                             the p-value for each metabolite will be calculated using 
                                                a t-test. The user can specificy log2FC and p-value thresholds 
                                                to establish significantly changed metabolites which will be used as 
                                                input for the pathway analysis (step 10)."),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("9. Identifier mapping")),
                                             h5("The metabolomics data is annotated with HMDB identifiers 
                                             , which can be transformed into ChEBI IDs using the ",
                                                tags$a(href = "https://bioconductor.org/packages/release/bioc/html/BridgeDbR.html", 
                                                       "BridgeDbR"), 
                                                       " (v2.6.0) package. This requires the .bridge file
                                                       which can be downloaded ",
                                                tags$a(href = "https://figshare.com/ndownloader/files/26001794", "here"),
                                                       " and then be placed in the ", em("data"), " folder.
                                                       Furthermore, in the applied procedure, one-to-many mappings 
                                                       are omitted by selecting the first hit only containing the prefix ‘CHEBI:’"),
                                             br(),
                                             hr(),
                                             
                                             
                                             h3(strong("10. Pathway analysis")),
                                             h5("The significantly changed metabolites from the 
                                             statistical analysis are retrieved as HMDB IDs and compared 
                                             against the content of each pathway in the WikiPathways SPARQL 
                                             endpoint using a semantic web query, filtering out pathways 
                                             from the Reactome database. Enrichment analysis of the 
                                             metabolomics data was conducted by over-representation 
                                             analysis (ORA) using a Fisher's 
                                             exact test by means of a hypergeometric density calculation.")
                                    ),
                                    
                                    #***************************************************#
                                    # Documentation of multi-omics visualization
                                    #***************************************************#
                                    tabPanel("Multi-omics Visualization",value = "doc_multi",
                                             h2(strong("Documentation of Multi-omics Visualization")),
                                             hr(),
                                             h3(strong("11. Pathway selection")),
                                             h5("In this step, the user can find interesting pathways by
                                                looking at the overlap of the significant pathways from 
                                                the transcriptomics and metabolomics analysis. The user 
                                                can select different significance thresholds for the metabolomics 
                                                and transcriptomics results, giving more weight to the 
                                                user-desired data type."),
                                             br(),
                                             hr(),
                                             h3(strong("12. Multi-omics visualization")),
                                             h5("This part requires Cytoscape to be installed and opened.
                                                The user-selected pathway will be loaded into Cytoscape 
                                                using the ",
                                                tags$a(href = "https://bioconductor.org/packages/release/bioc/html/RCy3.html", 
                                                       "RCy3"), 
                                                " and ",
                                                tags$a(href = "https://www.bioconductor.org/packages/release/bioc/html/rWikiPathways.html", 
                                                       "rWikiPathways"), 
                                                " (v1.16.0) packages and ",
                                                tags$a(href = "https://apps.cytoscape.org/apps/wikipathways", 
                                                       "WikiPathways app for Cytoscape"),
                                                " (v3.3.10). The gene expression 
                                                and metabolite data were visualized on the nodes of the 
                                                Cytoscape network as log2FC and p-values, based on their 
                                                Ensembl and ChEBI identifiers, respectively.")
                                             ),
                        )# tabsetpanel
               ) #tabPanel
               
               
               
    ) #navbar page
    
  )#fluidPage
)#tagList


