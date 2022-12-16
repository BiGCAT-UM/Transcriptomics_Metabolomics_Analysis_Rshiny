server = function(input, output,session) { 
  
  #to limit file size
  options(shiny.maxRequestSize = 125*1024^2)
  
  ###########Updating Side bar based on selected process##############
 
  #####Preprocess tab selected########
  output$trans_preprocess <- renderUI({
    tagList(
      conditionalPanel(condition = 'input.tabs=="Transcriptomics Analysis" && input.tabs_trans=="Preprocessing"',
                       
                       #Go forward
                       #******************************************************#
                       #   Information
                       #******************************************************#
                       
                       h3(strong("Pre-processing")),
                       h5("Sample filtering based on visit number and data type will be performed."),
                       br(),
                       h5("Samples with all zero values across all genes will be filtered."),
                       br(),
                       h5("Genes with all zero values across all samples will be filtered."),
                       tags$br(),
                       actionBttn(inputId ="filtering", label ="Apply", style = "jelly",
                                  btn_type = "button", type = "primary"),
                       tags$hr(),
                       h5("Average log CPM filtering will be performed to filter lowly expressed genes."),
                       textInput("threshold", "Threshold", 1),
                       verbatimTextOutput("value"),
                       actionBttn(inputId ="cpm_filtering", label ="Apply", style = "jelly",
                                  btn_type = "button", type = "primary"),
                       #Go forward
                       tags$br(),
                       tags$br(),
                       # Horizontal line ----
                       tags$hr(),
                       actionBttn(inputId ="preprocess_NEXT", label ="NEXT", style = "jelly",
                                  btn_type = "button", type = "primary")
          
      ))
  })#renderUI
  
  
  #####Data upload selected########
  output$trans_upload <- renderUI({
   
     tagList(
      conditionalPanel(condition = 'input.tabs=="Transcriptomics Analysis" && input.tabs_trans=="Data Upload"',
                       
                       fileInput("file1", "Choose a metadata file",
                                 multiple = FALSE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       # Input: Checkbox if file has header
                       checkboxInput("header1", "Header", TRUE),
                       # Input: Checkbox if file has row names
                       checkboxInput("rowNames1", "Rownames", TRUE),
                       # Input: Select separator ----
                       radioButtons("sep1", "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = ","),
                       # Horizontal line ----
                       tags$hr(),
                       
                       ###### -------------------COUNT DATA FILE READ --------####
                       # Input: Select a file ----
                       fileInput("file2", "Choose host transcriptomics count file",
                                 multiple = FALSE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       
                       # Input: Checkbox if file has header ----
                       checkboxInput("header2", "Header", TRUE),
                       # Input: Checkbox if file has row names
                       checkboxInput("rowNames2", "RowNames", TRUE),
                       # Input: Select separator ----
                       radioButtons("sep2", "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = "\t"),
                       
                       # Horizontal line ----
                       tags$hr(),
                       #Go forward
                       actionBttn(inputId ="upload_NEXT", label ="NEXT", style = "jelly",
                                 btn_type = "button", type = "primary")
                       
      )#conditionalPanel
      )#tagList
    
  }) #eof renderUI
  
  #####DEG analysis tab selected########
  output$trans_deg <- renderUI({
    tagList(
      conditionalPanel(condition = 'input.tabs=="Transcriptomics Analysis" && input.tabs_trans=="DEG Analysis"',
                      
                       h3(strong("Differential Gene Expression Analysis")),
                       h5("Statistical Analysis will be performed"),
                       br(),
                       
                       ########  Remove outliers #########
                       h4(strong("1. Remove samples")),
                       
                       awesomeCheckbox(inputId = "outlierCheckBox",
                                       label = "Keep all samples",
                                       value = TRUE,
                                       status = "danger"),
                       
                       uiOutput("outliersout"),
                       
                       h4(strong("2. Normalization & QC plots")),
                       
                       actionBttn(inputId ="outlierButton", label ="Apply", style = "jelly",
                                  btn_type = "button", type = "primary"),
                       hr(),
                       
                       ############################ DEG Analysis ###########################
                       h4(strong("3. DEG Analysis")),
                           # sliderInput(
                           #   inputId = "pthreshold",
                           #   label = "P threshold",
                           #   value = 0.05,
                           #   min = 0,
                           #   max = 0.5,
                           #   step = 0.01
                           # ),
                           
                           sliderInput(
                             inputId = "FCthreshold",
                             label = "FC threshold",
                             value = 0,
                             min = 0,
                             max = 5,
                             step = 0.1
                           ),
                       actionBttn(inputId ="DEGButton", label ="Apply", style = "jelly",
                                  btn_type = "button", type = "primary"),
                       hr(),
                       
                       h4(strong("4. Volcano Plots")),
                       actionBttn(inputId ="volcanoButton", label ="Apply", style = "jelly",
                                  btn_type = "button", type = "primary"),
                       
                       #Go forward
                       tags$br(),
                       tags$br(),
                       # Horizontal line ----
                       tags$hr(),
                       actionBttn(inputId ="deg_NEXT", label ="NEXT", style = "jelly",
                                  btn_type = "button", type = "primary")
                       
      ))
  })
  
  #####Identifier mapping tab selected########
  output$trans_mapping <- renderUI({
    tagList(
      conditionalPanel(condition = 'input.tabs=="Transcriptomics Analysis" && input.tabs_trans=="Identifier Mapping"',
                       
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
                       
      ))
  })
  
  #####Data upload selected########
  output$met_upload <- renderUI({
    
    tagList(
      conditionalPanel(condition = 'input.tabs=="Metabolomics Analysis" && input.tabs_trans=="Data Upload"',
                       
                       fileInput("file1", "Choose a metadata file",
                                 multiple = FALSE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       # Input: Checkbox if file has header
                       checkboxInput("header1", "Header", TRUE),
                       # Input: Checkbox if file has row names
                       checkboxInput("rowNames1", "Rownames", TRUE),
                       # Input: Select separator ----
                       radioButtons("sep1", "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = ","),
                      
                       # Horizontal line ----
                       tags$hr(),
                      
                       #Go forward
                       actionBttn(inputId ="metUpload_NEXT", label ="NEXT", style = "jelly",
                                  btn_type = "button", type = "primary")
                       
      )#conditionalPanel
    )#tagList
    
  })
  
  output$met_preprocess <- renderUI({
    
  })
  
  output$met_statistical<- renderUI({
    
  })
  
  
  ################# SAMPLE-GENE FILTERING BUTTON CLICK ################
 # FilterFunc<-reactive({sample_gene_filtering(input$metaData, input$countData)})
  
   data <- eventReactive(input$filtering, {
          
     sample_gene_filtering(metaData, countData )
                            
  })#eventReactive
   
  
   observeEvent(input$filtering, {
    
      showModal(modalDialog(
        title = "Message",
        paste0("Preprocessing started. Please wait it to be finished!"),
        easyClose = TRUE,
        footer = NULL
      ))
      
      #show preprocessed data 
      output$metaPreprocessed <- DT::renderDataTable(data()[[1]], server=TRUE)
      output$countPreprocessed <- DT::renderDataTable(data()[[2]], server=TRUE)
    
      showModal(modalDialog(
        title = "Process status",
        paste0("Samples and genes were filtered"),
        easyClose = TRUE,
        footer = NULL
      ))
      
  })#eof observeEvent
   
   
  ################ CPM FILTERING BUTTON CLICK #####################    
  observeEvent(input$cpm_filtering, {
    
    #Show loading message
    showModal(modalDialog(title = h4("CPM filtering is being applied....", 
                                     align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", 
                             align = "center"), easyClose = TRUE,
                          br()))
    
     output$plot <- renderPlot({
       cpm_filter( data()[[1]],  data()[[2]], input$threshold)
       
    
       
   })

  })
  
   
   
  #############################    DATA UPLOAD   ###########################
  observeEvent(input$upload_NEXT, {
    
    showModal(modalDialog(
      title = "Message",
      paste0("You can start analysis now!"),
      easyClose = TRUE,
      footer = NULL
    ))
    
     updateTabsetPanel(session, "tabs_trans",
                     selected = "Preprocessing"
                      )
    
  })#eof observeEvent
  
  output$fileContent1 <- DT::renderDataTable({
    req(input$file1)
    metaData <<- read.csv(input$file1$datapath,sep = input$sep1)
    return (metaData)
  }, server=TRUE)
  
  output$fileContent2<- DT::renderDataTable({
    req(input$file2)
    countData <<- read.csv(input$file2$datapath,sep = input$sep2)
    return (countData)
  }, server = TRUE)
  
 #######################    DATA UPLOAD ####################################
  
  
  
  ########################### DEG analysis  ###############################

  observeEvent(input$preprocess_NEXT, { 
    
    #Show loading message
    showModal(modalDialog(title = h4("You can apply DEG analysis now", 
                                     align = "center"), 
                          footer = NULL,
                          easyClose = TRUE,
                          br()))
    #show DEG analysis tab panel
     updateTabsetPanel(session, "tabs_trans",
                      selected = "DEG Analysis"
    )
    
  })#observeEvent
  
  
  #-----Select outlier samples for removal---
  output$outliersout <- renderUI({

    if (length(input$outlierCheckBox) > 0){
      if(input$outlierCheckBox == FALSE){
        #samples <- meta()[,1]
        pickerInput(inputId = "outliersPicker",
                    label = "Select samples to be removed",
                    choices = as.vector(colnames(data()[[2]])),
                    multiple = TRUE)

      }
    }
  })#renderUI


  
########## Normalize and QC plots and remove outliers  #######
  observeEvent(input$outlierButton, {
    

    ### if remove outlier not selected
    if (input$outlierCheckBox == TRUE){
       
       outliers <- NULL
          
       showModal(modalDialog(title = h4("Normalization and QC plots are
                                         being processing", 
                                        align = "center"), 
                             footer = NULL,
                             easyClose = TRUE,
                             br()))
       
        #the preprocessed data will be normalized
         normalize_QCplots(data()[[1]],  data()[[2]])
       
         showModal(modalDialog(title = h4("Normalization and QC plots done\n
                                          You can find all QC plots form 2-differential_gene_expression_analysis", 
                                           align = "center"), 
                                footer = NULL,
                                easyClose = TRUE,
                                br()))
          
    }

    ### if outlier removal selected
    if (input$outlierCheckBox == FALSE){
          
      if (length(input$outliersPicker) < 1){
            outliers <- NULL
          }
          if (length(input$outliersPicker) > 0){
            outliers <- input$outliersPicker
            
           #remove outliers
           data <-  removeOutliers( data()[[1]],  data()[[2]], outliers)
                            
            
           showModal(modalDialog(title = h4("Selected outlier(s) removed", 
                                             align = "center"), 
                                  footer = NULL,
                                  easyClose = TRUE,
                                  br()))
           
           showModal(modalDialog(title = h4("Normalization and QC plots are
                                         being processing", 
                                            align = "center"), 
                                 footer = NULL,
                                 easyClose = TRUE,
                                 br()))
           #normalize QC plots   
           normalize_QCplots(data()[[1]],  data()[[2]])
           
           showModal(modalDialog(title = h4("Normalization and QC plots done\n
                                             You can find all QC plots form 2-differential_gene_expression_analysis", 
                                             align = "center"), 
                                  footer = NULL,
                                  easyClose = TRUE,
                                  br()))
            
          }
    }
    
    observe({
      
      cat ("PCA plots will be shown\n")
      WORK_DIR <- getwd()

      path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/PCAanalysis__biopsylocation2.png")
      cat ("image path =",path,"\n")

       output$pcaPlotRaw <- renderImage({
        list(src = path, contentType = 'image/png',width = "500px", height = "500px",
             alt = "This is alternate text")

      })
       path2 <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/PCAanalysis__biopsylocation2.png")
       cat ("image path =",path,"\n")

       output$pcaPlotNorm <- renderImage({
         list(src = path2, contentType = 'image/png',width = "500px", height = "500px",
              alt = "This is alternate text")

       })

       
    })#observe

  })#observeEvent

  FC_threshold <- reactive({
    input$FCthreshold
  })

 
  ######################## DEG BUTTOON ########################
  
  summaryTable <- eventReactive(input$DEGButton, {
    
    DE_analysis  (data()[[1]],  data()[[2]], FC_threshold)
    
  })#eventReactive
  
  
  observeEvent(input$DEGButton, {
    
    showModal(modalDialog(
      title = "Message",
      paste0("DEG analysis started. Please wait it to be finished!"),
      easyClose = TRUE,
      footer = NULL
    ))
    
    #show preprocessed data 
    #output$outTable <- DT::renderDataTable(summaryTable(), server=TRUE)
    summaryTable()
  
    
    
    showModal(modalDialog(
      title = "Process status",
      paste0("Samples and genes were filtered"),
      easyClose = TRUE,
      footer = NULL
    ))
    
  })#eof observeEvent

  
  selectedValue <- reactiveVal()
  
  
   selectedRow <- eventReactive(input$compList,{
     print("compList_rows_selected")

   })

  
  ############################################################################
  
  
  observeEvent(input$volcanoButton, {
    
    showModal(modalDialog(title = h4("volcano plot started. Please wait...",
                                     align = "center"), 
                          footer = NULL,
                          easyClose = TRUE,
                          br())
    )
    
    #call the function for list showing
    df <- showFileList()
    
    # browser()
    
    output$compList   <- DT::renderDataTable (
      #df,selection = "single", 
      df, server = TRUE
    )
    
    #show selected file name
    output$selected <- renderText({
     selectedValue(" textOutput(, )")
      
    })
    
  
     
  showModal(modalDialog(title = h4("volcano plot finished\n
                                    You can find all plots under 2-differential_gene_expression_analysis",
                                     align = "center"), 
                          footer = NULL,
                          easyClose = TRUE,
                          br()))
    
  })#
  
  
  
}#eof server