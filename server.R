server = function(input, output,session) { 
  
  #to limit file size
  options(shiny.maxRequestSize = 125*1024^2)
  
  ################################################################
  
  # Transcriptomics
  
  # ################################################################
  # 
  hideTab("tabs_trans", target = "filtering_trans")
  hideTab("tabs_trans", target = "norm_trans")
  hideTab("tabs_trans", target = "deg_trans")
  hideTab("tabs_trans", target = "mapping_trans")
  hideTab("tabs_trans", target = "pathway_trans")
  hideTab("tabs_trans", target = "heatmap_trans")
  hideTab("tabs_trans", target = "network_trans")
  # 
  # ################################################################
  # 
  # # Metabolomics
  # 
  # ################################################################
  # 
  hideTab("tabs_mets", target = "filtering_mets")
  hideTab("tabs_mets", target = "stat_mets")
  hideTab("tabs_mets", target = "pathway_mets")
  hideTab("tabs_mets", target = "mapping_mets")

  #**************************************************************************************************************#
  #                      Transcriptomics Data Operations
  #**************************************************************************************************************#
  
  #***************************************************#
  # Data Upload
  #***************************************************#
  
  # Load meta data
  metaData1 <- reactive({
    req(input$file1)
    metaData1 <- read.csv(input$file1$datapath,sep = input$sep1)
    return (metaData1)
  })
  
  # Return table with meta data
  output$fileContent1 <- DT::renderDataTable({
    metaData1()
  }, server=TRUE, options = list(pageLength = 5), rownames= FALSE)
  
  # Return header for meta data
  observeEvent((length(metaData1())> 0),{
    output$metaText <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })
  
  
  # Load transcriptomics count data
  countData1 <- reactive({
    req(input$file2)
    countData1 <- read.csv(input$file2$datapath,sep = input$sep2)
    return (countData1)
  })
  
  # Return table with transcriptomics data
  output$fileContent2 <- DT::renderDataTable({
    countData1()
  }, server = TRUE, options = list(pageLength = 5), rownames= TRUE)
  
  # Return header for transcriptomics data
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
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully selected! You can now start with the pre-processing!",
      type = "success")
    
    # Go to next tab
    updateTabsetPanel(session, "tabs_trans",
                      selected = "filtering_trans")
    
    # Show next tab
    showTab("tabs_trans", target = "filtering_trans")
    
  })
  
  
  #***************************************************#
  # Data filtering
  #***************************************************#
  
  # Perform Sample/Gene filtering
  data <- eventReactive(input$filtering, {
    data <- sample_gene_filtering(metaData1(), countData1())
    return(data)
  })
  
  # Return table with filtered meta data
  output$metaPreprocessed <- DT::renderDataTable(data()[[1]], server=TRUE,
                                                 options = list(pageLength = 5))
  
  
  # Return header for meta data
  observeEvent((length(data())> 0),{
    output$metaText1 <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })
  
  
  # Return table for transcriptomics data
  output$countPreprocessed <- DT::renderDataTable(data()[[2]], server=TRUE,
                                                  options = list(pageLength = 5))
  
  # Return header for transcriptomics data
  observeEvent((length(data())> 0),{
    output$countText1 <- renderUI({
      tagList(
        br(),
        hr(),
        h3(strong("Transcriptomics data"))
      )
    })
  })
  
  # Send success message
  observeEvent(input$filtering, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Samples and genes were successfully filtered!",
      type = "success")
    
  })
  
  # Render no plot
  output$plot <- renderPlot({
    NULL     
  })
  
  
  # Perform logCPM filtering
  data_filtered <- eventReactive(input$filtering, {
    data_filtered <- cpm_filter_output(data()[[1]],data()[[2]], input$threshold) 
    return(data_filtered)
  })
  
  # Header for Log CPM filtering   
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
    
    # Plot log CPM giltering
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
  observeEvent(if (length(data_filtered())> 0){input$preprocess_NEXT}, {
    
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
      normalize_QCplots(data_filtered()[[1]], data_filtered()[[2]])
      
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
        data1 <-  removeOutliers(data_filtered()[[1]], data_filtered()[[2]], outliers)
        
        
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
        
        list(src = path, contentType = 'image/png',width = "800px", height = "auto",
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
  
  P_threshold <- reactive({
    input$pthreshold
  })
  
  # Read top table
  topTable <- eventReactive(input$DEGButton, {
    req(input$AdjOrRaw)
    topTable <- list()
    WORK_DIR <- getwd()
    
    if (input$AdjOrRaw == "Raw") {
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab"))
      topTable[[1]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab"))
      topTable[[2]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab"))
      topTable[[3]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      
      temp<- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab"))
      topTable[[4]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
    }
    if (input$AdjOrRaw == "Adjusted") {
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab"))
      topTable[[1]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab"))
      topTable[[2]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      temp <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab"))
      topTable[[3]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      
      temp<- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab"))
      topTable[[4]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
    }
    

    
    names(topTable) <- c("Ileum: CD vs non-IBD",
                         "Rectum: CD vs non-IBD",
                         "Ileum: UC vs non-IBD",
                         "Rectum: UC vs non-IBD")
    
    return(topTable)
  })
  
  output$topTable <- DT::renderDataTable({
    req(input$Comparison)
    output <- topTable()[[input$Comparison]]
    output <- arrange(output, pvalue)
    return(output)
  }, server=TRUE,
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
  
  #**************************************************************************************************************#
  #                      Metabolomics Data Operations
  #**************************************************************************************************************#
  
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
  
  
  
  #output of metabolite intensity data
  observeEvent( input$metDownload,{
    
    #download and read metabolomics peak intensity data
    if(file.exists("data/metabolomics.csv.gz"))
      {print("Metabolomics zipped data already downloaded")}
    else{
      fileUrl <- "https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz?accessType=DOWNLOAD"
      download(fileUrl, "data/metabolomics.csv.gz", mode = "wb")
    }
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "You now have metabolomics count data! You can now start analysis!",
      type = "success")
   
  })#observeEvent

  #go to next step
  observeEvent(if ((length(mbxMeta())> 0)){input$metUpload_NEXT}, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully uploaded! You can now start with the pre-processing!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_mets",
                      selected = "filtering_mets")
    
    showTab("tabs_mets", target = "filtering_mets")
    
  })#eof observeEvent
 
  #***************************************************#
  # Data preprocessing
  #***************************************************#
  #*
  # Metabolomics Sample/Gene filtering
  data2 <- eventReactive(input$metsFiltering, {
    sendSweetAlert(
      session = session,
      title = "Message",
      text = "Filtering process started! It might take time please be patient.",
      type = "info"
      )
    preprocessMets(mbxMeta())
    
  })#eventReactive


  # Filtered meta data
  observeEvent((length(data2())> 0),{
    output$metaPreprocessText <- renderUI({
      tagList(
        h3(strong("Meta data"))
      )
    })
  })

  output$mbxMetaPreprocessed <- DT::renderDataTable(data2()[[1]], server=TRUE,
                                                 options = list(pageLength = 5))

  # Filtered  data
  observeEvent((length(data2())> 0),{
    output$mbxCountText <- renderUI({
      tagList(
        br(),
        hr(),
        h3(strong("Metabolomics data"))
      )
    })
  })
  
  output$mbxCountPreprocessed <- DT::renderDataTable(data2()[[2]], server=TRUE,
                                                  options = list(pageLength = 5))

  observeEvent(input$metsFiltering, {

    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Samples and genes were successfully filtered!",
      type = "success")

  })#eof observeEvent
  
  
  
  # Go the next step
  observeEvent(if ((length(data2()[[1]])> 0) && (length(data2()[[2]])> 0)){input$metPre_NEXT}, {

    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully preprocessed! You can now start with statistical analysis!",
      type = "success")

    updateTabsetPanel(session, "tabs_mets",
                      selected = "stat_mets")

    showTab("tabs_mets", target = "stat_mets")

  })#eof observeEvent
  
  
  
  
}#eof server