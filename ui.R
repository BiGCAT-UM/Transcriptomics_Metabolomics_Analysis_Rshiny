ui <- tagList(
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
                                    # Preprocessing
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
                                               h4(strong("Average log CPM filtering")),
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
                                                 label = "p-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
                                               ),
                                               br(),
                                               numericInput(
                                                 inputId = "q_threshold_pathway",
                                                 label = "q-value threshold",
                                                 value = 0.02,
                                                 min = 0,
                                                 max = 1,
                                                 step = 0.01
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
                                                  the analysis"),
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
                                                           choices = c("ileum", "rectum"),
                                                           selected = "ileum"),
                                               imageOutput("NetworkPlot",
                                                           width = "700px",
                                                           height = "auto")
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
                                               hr(),
                                               h5("The peak intensity value data will be uploaded,"),
                                               
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
                                                            selected = "\t"),
                                               
                                               # Horizontal line
                                               tags$hr(),
                                               
                                               # Horizontal line ----
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
                                               h5("1. Samples based on visit number and data type will be filtered."),
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
                                               h4("Please select a normalization method to be applied:"),
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
                                    
                                    
                                    #Statistical Analysis
                                    #***************************************************#
                                    tabPanel("Statistical Analysis", value= "stat_mets",

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
                                          value = 1.5,
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
                                               h5("To visualize the multi-omics data in CytoScape, data need to be prepared 
                                               for importing according to the unified database identifiers available.
                                               Therefor, metabolite data is mapped to ChEBI IDS."),
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
                                               DT::dataTableOutput("mappingTable")
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
                                          
                                          # Horizontal line ----
                                          tags$hr(),
                                          actionBttn(inputId ="pathwayMet_NEXT", label ="Next", style = "jelly",
                                                     btn_type = "button", type = "primary", color = "danger",
                                                     icon = icon("arrow-right"))
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
               
               
               # # Multi-omics
               # 
               # ################################################################
               tabPanel("Multi-omics Visualization",
                        value = "multiomics",
                        tabsetPanel(id="tabs_multi",

                                    #***************************************************#
                                    # Pathway Selection
                                    #***************************************************#
                                    tabPanel("Pathway Selection",value = "selection"),

                                    #***************************************************#
                                    # Visualization
                                    #***************************************************#
                                    tabPanel("Visualization", value = "visualization")

                        )#tabsetPanel
               )#tabPanel
               
               
               
    ) #navbar page
    
  )#fluidPage
)#tagList


