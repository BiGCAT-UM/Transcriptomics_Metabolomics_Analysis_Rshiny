
  ui = tags$body(class="skin-blue sidebar-mini control-sidebar-open",
                 dashboardPage(
                   options = list(sidebarExpandOnHover = TRUE),
                   header = dashboardHeader(title = "Functional Analysis of transcriptomics and metabolomics",titleWidth = 450),
                   
                   sidebar = dashboardSidebar(minified = F, collapsed = F,
                                              h4("Selected Processes"),
                                              uiOutput("trans_upload"),uiOutput("trans_preprocess"),uiOutput("trans_deg"),
                                              uiOutput("trans_mapping"),uiOutput("trans_pathway"),uiOutput("trans_heatmap"),uiOutput("trans_network")
                                              
                   ),#eof side bar
                   
                   
                   body = dashboardBody(
                     
                     tabsetPanel(id = "tabs",
                                 tabPanel("Transcriptomics Analysis", 
                                          tabsetPanel(id="tabs_trans",
                                                      tabPanel("Data Upload",
                                                               mainPanel(
                                                                 # Output: Data file-1 ----
                                                                 DT::dataTableOutput("fileContent1"),
                                                                 # Output: Data file-2 ----
                                                                 DT::dataTableOutput("fileContent2")
                                                               )
                                                              ),
                                                      
                                                      tabPanel("Preprocessing",
                                                            
                                                               
                                                               mainPanel(
                                                                 # Output: Data file-1 ----
                                                                 DT::dataTableOutput("metaPreprocessed"),
                                                                 br(),
                                                                 # Output: Data file-2 ----
                                                                 DT::dataTableOutput("countPreprocessed"),
                                                                 br(),
                                                                 plotOutput("plot")
                                                               )
                                                        ),
                                                      
                                                      
                                                      tabPanel("DEG Analysis" ,
                                                         
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
                                                            
                                                          )#fluidPage
                                                          
                                                         
                                                         ),#mainPanel
                                                         
                                                         height = "1000px",
                                                         width = "1000px"
                                                      ),#tabPanel
                                                      
                                                      tabPanel("Identifier Mapping",
                                                               mainPanel(
                                                                 # Output: -
                                                                 
                                                               )
                                                               
                                                               
                                                               
                                                      ),
                                                      tabPanel("Pathway Analysis"),
                                                      
                                                      tabPanel("Create Heatmap"),
                                                      
                                                      tabPanel("Network Analysis")
                                                      
                                          ),
                                          
                                         
                                          
                                 ),#eof tabPanel
                                 
                                 tabPanel("Metabolomics Analysis",
                                          tabsetPanel(id="tabs_mets",
                                                      tabPanel("Preprocessing",value = "pre_mets"),
                                                      tabPanel("Statistical Analysis", value= "stat_mets"),
                                                      tabPanel("Pathway Analysis", value = "pathway_mets"),
                                                      tabPanel("Identifier Mapping", value = "mapping_mets")
                                                      
                                          )
                                          
                                 ),#eof tabPanel
                                 
                                 tabPanel("Multi-omics Visualization",
                                          tabsetPanel(id="tabs_multi",
                                                      tabPanel("Pathway Selection",value = "selection"),
                                                      tabPanel("Visualization", value = "visualization")
                                                      
                                          )
                                          
                                 )#eof tabPanel
                                 
                     ),#eof tabsetPanel   
                     
                   
                     
                     
                   ),#eof dashboard
                   
                   
                   title = "DashboardPage"
                 )#eof dashboard
                 
  )
  

