server = function(input, output,session) { 
  
  #to limit file size
  options(shiny.maxRequestSize = 125*1024^2)
  
  ################################################################
  
  # Transcriptomics
  
  ################################################################
  
  hideTab("tabs_trans", target = "filtering_trans")
  hideTab("tabs_trans", target = "norm_trans")
  hideTab("tabs_trans", target = "deg_trans")
  hideTab("tabs_trans", target = "mapping_trans")
  hideTab("tabs_trans", target = "pathway_trans")
  hideTab("tabs_trans", target = "heatmap_trans")
  hideTab("tabs_trans", target = "network_trans")
  
  ################################################################
  
  # Metabolomics
  
  ################################################################
  
  hideTab("tabs_mets", target = "filtering_mets")
  hideTab("tabs_mets", target = "stat_mets")
  hideTab("tabs_mets", target = "pathway_mets")
  hideTab("tabs_mets", target = "mapping_mets")
  
  
  #***************************************************#
  # Data Upload
  #***************************************************#
  
  # Output of meta data
  metaData1 <- reactive({
    req(input$file1)
    metaData1 <- read.csv(input$file1$datapath,sep = input$sep1)
    return (metaData1)
  })
  
  output$fileContent1 <- DT::renderDataTable({
    metaData1()
  }, server=TRUE, options = list(pageLength = 5), rownames= FALSE)
  
  
  observeEvent((length(metaData1())> 0),{
    output$metaText <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })
  
  
  # Output of transcriptomics data
  countData1 <- reactive({
    req(input$file2)
    countData1 <- read.csv(input$file2$datapath,sep = input$sep2)
    return (countData1)
  })
  
  output$fileContent2 <- DT::renderDataTable({
    countData1()
  }, server = TRUE, options = list(pageLength = 5), rownames= FALSE)
  
  
  observeEvent((length(countData1())> 0),{
    output$countText <- renderUI({
      tagList(
        br(),
        hr(),
        h3(strong("Transcriptomics data"))
      )
    })
  })
  
  # Go the next step
  observeEvent(if ((length(countData1())> 0) && (length(metaData1())> 0)){input$upload_NEXT}, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully selected! You can now start with the pre-processing!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "filtering_trans")
    
    showTab("tabs_trans", target = "filtering_trans")
    
  })#eof observeEvent
  
  
  #***************************************************#
  # Data filtering
  #***************************************************#
  
  
  # Sample/Gene filtering
  data <- eventReactive(input$filtering, {
    sample_gene_filtering(metaData1(), countData1())
  })#eventReactive
  
  # Filtered meta data
  observeEvent((length(data())> 0),{
    output$metaText1 <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })
  
  output$metaPreprocessed <- DT::renderDataTable(data()[[1]], server=TRUE,
                                                 options = list(pageLength = 5))
  
  # Filtered transcriptomics data
  observeEvent((length(data())> 0),{
    output$countText1 <- renderUI({
      tagList(
        br(),
        hr(),
        h3(strong("Transcriptomics data"))
      )
    })
  })
  
  output$countPreprocessed <- DT::renderDataTable(data()[[2]], server=TRUE,
                                                  options = list(pageLength = 5))
  
  observeEvent(input$filtering, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Samples and genes were successfully filtered!",
      type = "success")
    
  })#eof observeEvent
  
  output$plot <- renderPlot({
    NULL     
  })
  
  # Log CPM filtering   
  observeEvent(input$cpm_filtering, {
    observeEvent(if (length(data())> 0){input$cpm_filtering},{
      output$histText <- renderUI({
        tagList(
          br(),
          hr(),
          h3(strong("Histogram of Mean Expression Values"))
        )
      })
    })
    
    output$plot <- renderPlot({
      cpm_filter(data()[[1]],data()[[2]], input$threshold)       
    })
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Log CPM filtering was successfully applied",
      type = "success")
    
  })
  
  # Go the next step
  observeEvent(if (length(data())> 0){input$preprocess_NEXT}, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully filtered! You can now start with the normalization!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "norm_trans")
    
    showTab("tabs_trans", target = "norm_trans")
    
  })
  
  #***************************************************#
  # Normalization and QC
  #***************************************************#
  
  
  #====================================================#
  # Select outlier(s)
  #====================================================#
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
  })
  
  #====================================================#
  # Normalize data
  #====================================================#
  observeEvent(input$outlierButton, {
    
    ### if remove outlier not selected
    if (input$outlierCheckBox == TRUE){
      
      outliers <- NULL
      
      showModal(modalDialog(title = h4(strong("Normalization and Quality Control"),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      #the preprocessed data will be normalized
      normalize_QCplots(data()[[1]], data()[[2]])
      
      removeModal()
      
      sendSweetAlert(
        session = session,
        title = "Success!",
        text = "Normalization and QC plots done. 
        You can find all QC plots in 2-differential_gene_expression_analysis",
        type = "success")
      
    }
    
    ### if outlier removal selected
    if (input$outlierCheckBox == FALSE){
      
      if (length(input$outliersPicker) < 1){
        outliers <- NULL
      }
      if (length(input$outliersPicker) > 0){
        outliers <- input$outliersPicker
        
        #remove outliers
        data1 <-  removeOutliers(data()[[1]], data()[[2]], outliers)
        
        
        showModal(modalDialog(title = h4(strong("Normalization and Quality Control"),
                                         align = "center"), 
                              footer = NULL,
                              h5("This might take a while. Please be patient.", 
                                 align = "center")))
        
        #normalize QC plots   
        normalize_QCplots(data1[[1]], data1[[2]])
        
        removeModal()
        
        sendSweetAlert(
          session = session,
          title = "Success!",
          text = "Normalization and QC plots done. 
        You can find all QC plots in 2-differential_gene_expression_analysis",
          type = "success")
        
      }
    }
    
    
    #====================================================#
    # Make PCA plots
    #====================================================#
    observe({
      
      cat ("PCA plots will be shown\n")
      WORK_DIR <- getwd()
      
      
      output$QCplot <- renderImage({
        req(input$whichQCplot)
        if (input$whichQCplot == "PCA (Normalized)"){
          path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/PCAanalysis__biopsylocation2.png")
          cat ("image path =",path,"\n")
        }
        if (input$whichQCplot == "PCA (Raw)"){
          path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/PCAanalysis__biopsylocation2.png")
          cat ("image path =",path,"\n")
        }
        if (input$whichQCplot == "Boxplot (Normalized)"){
          path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/Boxplot__biopsylocation.png")
          cat ("image path =",path,"\n")
        }
        if (input$whichQCplot == "Boxplot (Raw)"){
          path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/Boxplot__biopsylocation.png")
          cat ("image path =",path,"\n")
        }
        
        list(src = path, contentType = 'image/png',width = "700px", height = "auto",
             alt = "This is alternate text")
        
      }, deleteFile=FALSE)
      
      
    })#observe
    
  })#observeEvent
  
  # Go the next step
  observeEvent(input$norm_NEXT, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully normalized! You can now start with the DEG analysis!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "deg_trans")
    
    showTab("tabs_trans", target = "deg_trans")
    
  })
  
  #***************************************************#
  # DEG analysis
  #***************************************************#
  
  FC_threshold <- reactive({
    input$FCthreshold
  })
  

  topTable <- eventReactive(input$DEGButton, {
    WORK_DIR <- getwd()
    topTable <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab"))
    return(topTable)
  })
  
  output$topTable <- DT::renderDataTable(topTable(), server=TRUE,
                                                 options = list(pageLength = 5))
  
  observeEvent(input$DEGButton, {
    
    showModal(modalDialog(title = h4(strong("Statistical Analysis"),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    
    
    removeModal()
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "DEG analysis is finished!",
      type = "success")
    
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
  
  #***************************************************#
  #                      Metabolomics Data Operations
  #***************************************************#
  
  #***************************************************#
  # Data Upload
  #***************************************************#
  
  # Output of meta data
  mbxMeta <- reactive({
    req(input$metaFile)
    mbxMeta <- read.csv(input$metaFile$datapath,sep = input$sepMet)
    return (mbxMeta)
  })
  
  output$mbxMetaFile <- DT::renderDataTable({
    mbxMeta()
  }, server=TRUE, options = list(pageLength = 5), rownames= FALSE)
  
  
  observeEvent((length(mbxMeta())> 0),{
    output$mbxMetaText <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })
  
  
  
  
  
  
}#eof server