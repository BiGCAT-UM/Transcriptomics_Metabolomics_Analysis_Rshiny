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
                                               selectInput(inputId = "whichQCplot",
                                                           label = NULL,
                                                           choices = c("PCA (Normalized)",
                                                                       "PCA (Raw)",
                                                                       "Boxplot (Normalized)",
                                                                       "Boxplot (Raw)"),
                                                           selected = "PCA (Normalized)"),
                                               
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
                                                            choices = c("Adjusted",
                                                                        "Raw"),
                                                            selected = "Adjusted"),
                                               br(),
                                               numericInput(
                                                 inputId = "pthreshold",
                                                 label = "P-value threshold",
                                                 value = 0.05
                                               ),
                                               br(),
                                               
                                               numericInput(
                                                 inputId = "FCthreshold",
                                                 label = "FC threshold",
                                                 value = 0
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
                                               
                                               selectInput(inputId = "Comparison",
                                                           label = NULL,
                                                           choices = c("Ileum: CD vs non-IBD",
                                                                       "Rectum: CD vs non-IBD",
                                                                       "Ileum: UC vs non-IBD",
                                                                       "Rectum: UC vs non-IBD"),
                                                           selected = "Ileum: CD vs non-IBD"),
                                               
                                               #output from DE analysis
                                               DT::dataTableOutput("topTable")%>% 
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
                                               h5("The peak intensity value data will be downloaded,"),
                                               h5("If you already downloaded it, the count data will be uploaded"),
                                              
                                               #Metabolomics data download button
                                               actionBttn(inputId ="metDownload", 
                                                          label ="Download Data", 
                                                          style = "jelly",
                                                          btn_type = "button", 
                                                          type = "primary"),
                                               
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
                                                 withSpinner(color="#0dc5c1")
                                            
                                               
                                             )#mainPanel
                                             
                                    ),#tabPanel
                                    
                                    
                                    #***************************************************#
                                    # Preprocessing
                                    #***************************************************#
                                    tabPanel("Filtering",value = "filtering_mets",
                                             
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               # Title + description
                                               h3(strong("Pre-processing")),
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
                                               uiOutput("metaPreprocessText"),
                                               DT::dataTableOutput("mbxMetaPreprocessed")%>%
                                                 withSpinner(color="#0dc5c1"),
                                               br(),
                                               # Output: Data file-2 ----
                                               uiOutput("mbxCountText"),
                                               DT::dataTableOutput("mbxCountPreprocessed")%>%
                                                 withSpinner(color="#0dc5c1"),
                                               
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
                                               selectInput(inputId = "whichQCplot_met",
                                                           label = NULL,
                                                           choices = c("PCA (Normalized)",
                                                                       "PCA (Raw)",
                                                                       "Boxplot (Normalized)",
                                                                       "Boxplot (Raw)"),
                                                           selected = "PCA (Normalized)"),
                                               
                                               imageOutput("QCplot_met",
                                                           width = "700px",
                                                           height = "auto"),
                                               
                                               
                                               
                                               
                                               
                                             )#End of mainPanel
                                             
                                    ),
                                    
                                    
                                    
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
                        
               ),#tabPanel
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


