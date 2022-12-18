ui <- tagList(
  fluidPage(
    
    # Allow arning/information messages
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
                                               
                                               # Upload meta data
                                               fileInput(inputId = "file1", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv")),
                                               
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
                                               
                                               # Upload file
                                               fileInput(inputId = "file2", 
                                                         label = NULL,
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv")),
                                               
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
                                               h3(strong("Pre-processing")),
                                               hr(),
                                               h4(strong("Sample and gene filtering")),
                                               h5("1. Samples with all zero values across all genes will be filtered."),
                                               h5("2. Genes with all zero values across all samples will be filtered."),
                                               br(),
                                               
                                               # Apply filtering
                                               actionBttn(inputId ="filtering", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               tags$hr(),
                                               
                                               # Title + description
                                               h4(strong("Average log CPM filtering")),
                                               h5("This filtering step will be performed to remove lowly expressed genes."),
                                               
                                               # log CPM threshold
                                               numericInput(inputId = "threshold", 
                                                            label = "Log CPM threshold", 
                                                            value = 1),
                                               
                                               actionBttn(inputId ="cpm_filtering", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary", color = "primary"),
                                               
                                               #Go forward
                                               tags$br(),
                                               tags$br(),
                                               tags$hr(),
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
                                               # Output: Data file-1 ----
                                               uiOutput("metaText1"),
                                               DT::dataTableOutput("metaPreprocessed")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                               br(),
                                               # Output: Data file-2 ----
                                               uiOutput("countText1"),
                                               DT::dataTableOutput("countPreprocessed")%>% 
                                                 withSpinner(color="#0dc5c1"),
                                               br(),
                                               uiOutput("histText"),
                                               plotOutput("plot")%>% 
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
                                               imageOutput("pcaPlotNorm",
                                                           width = "1000px",
                                                           height = "800px"),
                                               
                                               
                                                 
                                               
                                               
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
                                                            choices = c("Adjusted",
                                                                        "Raw"),
                                                            selected = "Adjusted"),
                                               br(),
                                               sliderInput(
                                                 inputId = "pthreshold",
                                                 label = "P-value threshold",
                                                 value = 0.05,
                                                 min = 0,
                                                 max = 0.5,
                                                 step = 0.01
                                               ),
                                               br(),
                                               
                                               sliderInput(
                                                 inputId = "FCthreshold",
                                                 label = "log2FC threshold",
                                                 value = 0,
                                                 min = -5,
                                                 max = 5,
                                                 step = 0.5
                                               ),
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
                                               
                                               fluidPage(
                                                 
                                                 #output from DE analysis
                                                 tableOutput("outTable"),
                                                 
                                                 #outputs from volcano plots
                                                 fluidRow(
                                                   column(3,DTOutput("compList")),
                                                   br(),
                                                   verbatimTextOutput('selected'),
                                                   column(9,imageOutput("volcanoPlot"))
                                                 )
                                                 
                                               )#End of fluidPage
                                               
                                             ),#End of mainPanel
                                             
                                             height = "1000px",
                                             width = "1000px"
                                             
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
                                               h5("HGNC gene symbols will be transformed into ENTREZ IDs"),
                                               br(),
                                               
                                               actionBttn(inputId ="mapping", label ="Apply", style = "jelly",
                                                          btn_type = "button", type = "primary"),
                                               #Go forward
                                               tags$br(),
                                               tags$br(),
                                               # Horizontal line ----
                                               tags$hr(),
                                               actionBttn(inputId ="mapping_NEXT", label ="NEXT", style = "jelly",
                                                          btn_type = "button", type = "primary")
                                               
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               # Output: -
                                               
                                             )
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Pathway Analysis
                                    #***************************************************#
                                    
                                    tabPanel("Pathway Analysis", value = "pathway_trans"),
                                    
                                    #***************************************************#
                                    # Heatmap
                                    #***************************************************#
                                    
                                    tabPanel("Create Heatmap", value = "heatmap_trans"),
                                    
                                    #***************************************************#
                                    # Network analysis
                                    #***************************************************#
                                    
                                    tabPanel("Network Analysis", value = "network_trans")
                                    
                        ) # End of tabset pabel
                        
               ),
               
               ################################################################
               
               # Metabolomics
               
               ################################################################
               tabPanel("Metabolomics Analysis", 
                        value = "metabolomics",
                        
                        tabsetPanel(id="tabs_mets",
                                    
                                    #***************************************************#
                                    # Preprocessing
                                    #***************************************************#
                                    tabPanel("Preprocessing",value = "pre_mets"),
                                    
                                    #***************************************************#
                                    #Statistical Analysis
                                    #***************************************************#
                                    tabPanel("Statistical Analysis", value= "stat_mets"),
                                    
                                    #***************************************************#
                                    # Pathway Analysis
                                    #***************************************************#
                                    tabPanel("Pathway Analysis", value = "pathway_mets"),
                                    
                                    #***************************************************#
                                    # Identifier mapping
                                    #***************************************************#
                                    tabPanel("Identifier Mapping", value = "mapping_mets")
                                    
                        )
                        
               ),
               ################################################################
               
               # Multi-omics
               
               ################################################################
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
                                    
                        )
               )
               
               
               
    ) # End of navbar page
    
  )
)


