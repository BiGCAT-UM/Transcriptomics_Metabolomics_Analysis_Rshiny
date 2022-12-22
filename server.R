server = function(input, output,session) { 
  
  #to limit file size
  options(shiny.maxRequestSize = 500*1024^2)
  
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
  #hideTab("tabs_mets", target = "norm_mets")
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
  
  # Perform Sample/Gene filtering and logCPM filtering
  data <- eventReactive(input$cpm_filtering, {
    data <- sample_gene_filtering(metaData1(), countData1())
    return(data)
  })
  
  data_filtered <- reactive({
    req(data())
    data_filtered <- cpm_filter_output(data()[[1]],data()[[2]], input$threshold) 
    return(data_filtered)
  })
  
  output$plot_prefiltering <- renderPlot(NULL)
  output$plot_postfiltering <- renderPlot(NULL)
  
  # Header for Log CPM filtering   
  observeEvent(input$cpm_filtering, {
    
    observeEvent(if (length(data())> 0){input$cpm_filtering},{
      output$histText_pre <- renderUI({
        tagList(
          br(),
          hr(),
          h3(strong("Histogram of Mean Expression Values (Pre-filtering)"))
        )
      })
    })
    
    # Plot log CPM filtering
    output$plot_prefiltering <- renderPlot({
      cpm_filter(data()[[1]],data()[[2]], input$threshold)       
    })
    
    # Header for Log CPM filtering   
    observeEvent(if (length(data())> 0){input$cpm_filtering},{
      output$histText_post <- renderUI({
        tagList(
          br(),
          hr(),
          h3(strong("Histogram of Mean Expression Values (Post-filtering)"))
        )
      })
    })
    
    # Plot log CPM giltering
    output$plot_postfiltering <- renderPlot({
      cpm_filter(data_filtered()[[1]],data_filtered()[[2]], NULL)       
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
  
  # FC threshold
  FC_threshold <- reactive({
    input$FCthreshold
  })
  
  # P value threshold
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

  # Return top table
  output$topTable <- DT::renderDataTable({
    req(input$Comparison)
    output <- topTable()[[input$Comparison]]
    output <- arrange(output, pvalue)
    output <- output[,c(1,2,3,4,7,8)]
    colnames(output) <- c("Gene Name", "Avg Expr", "log2 FC", "FC", "p-value", "adj. p-value")
    return(output)
  }, server=TRUE,
  options = list(pageLength = 5), rownames= FALSE)
  
  output$VolcanoPlot <- renderPlot(NULL)
  
  observeEvent(input$DEGButton, {
    
    # Make plot
    output$VolcanoPlot <- renderPlot({
      req(input$AdjOrRaw)
      req(input$Comparison)
      topTable <- topTable()[[input$Comparison]]
      if (input$AdjOrRaw == "Adjusted") {
        p<-EnhancedVolcano(topTable, lab = topTable$X, labSize = 3, title = NULL, 
                           x = 'log2FoldChange',
                           y = 'pvalue', 
                           pCutoff = P_threshold(), 
                           FCcutoff = log2(FC_threshold()))
      }
      if (input$AdjOrRaw == "Raw") {
        p<-EnhancedVolcano(topTable, lab = topTable$X, labSize = 3, title = NULL, 
                           x = 'log2FoldChange',
                           y = 'padj', 
                           pCutoff = P_threshold(), 
                           FCcutoff = log2(FC_threshold()))
      }
      return(p)
    })
    
    # Loading message
    showModal(modalDialog(title = h4(strong("Statistical Analysis"),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    
    
    removeModal()
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "DEG analysis is finished!",
      type = "success")
    
  })#eof observeEvent
  

  
  
  
  
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
  
  # Output of metabolomics original data
  mbxData <- reactive({
    req(input$metsData)
    mbxData <- read.csv(input$metsData$datapath,sep = input$sepMet2)
    return (mbxData)
  })
  
  output$mbxDataFile <- DT::renderDataTable({
    mbxData()
  }, server=TRUE, options = list(pageLength = 5), rownames= FALSE)
  
  
  observeEvent((length(mbxData())> 0),{
    output$mbxDataFileText <- renderUI({
      tagList(
        h3(strong("Metabolomics data"))
      )
    })
  })
   ##################################################################

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
  # Data filtering
  #***************************************************#
  
  # Metabolomics Sample/Gene filtering
  data2 <- eventReactive(input$metsFiltering, {
    sendSweetAlert(
      session = session,
      title = "Message",
      text = "Filtering process started! It might take time please be patient.",
      type = "info"
      )
    data2 <- filteringMets(mbxMeta(),mbxData())
  return (data2)
  })#eventReactive


  # Header for filtered meta data
  observeEvent((length(data2())> 0),{
    output$CDpreprocessText <- renderUI({
      tagList(
        h3(strong("CD preprocessed metabolomics data"))
      )
    })
  })

  #Table for filtered meta-data
  output$mbxCDPreprocessed <- DT::renderDataTable(data2()[[1]], server=TRUE,
                                                 options = list(pageLength = 5))

  # Header for filtered metabolomics data
  observeEvent((length(data2())> 0),{
    output$UCpreprocessText <- renderUI({
      tagList(
        h3(strong("UC preprocessed metabolomics data"))
      )
    })
  })
  
  #Table for filtered metabolomics data
  output$mbxUCPreprocessed <- DT::renderDataTable(data2()[[2]], server=TRUE,
                                                  options = list(pageLength = 5))

  observeEvent(input$metsFiltering, {

    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully filtered!",
      type = "success")

  })#eof observeEvent
  
  # Go the next step
  observeEvent(if ((length(data2()[[1]])> 0) && (length(data2()[[2]])> 0)){input$metPre_NEXT}, {

    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully filtered! You can now start with normalization process!",
      type = "success")

    updateTabsetPanel(session, "tabs_mets",
                      selected = "norm_mets")

    showTab("tabs_mets", target = "norm_mets")

  })#eof observeEvent
  
  #####################################################################################    
  
  #***************************************************#
  # Data normalization
  #***************************************************#
 
  selectedMethod <- reactive({
    input$whichNormMethod
  })
  
  observeEvent( input$normButton,{
    
    #normalization process is done with the function below
    normalizeMets(selectedMethod())
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Normalization and QC plots done. 
        You can find all QC plots in 7-metabolite",
      type = "success")
    
    cat ("Histograms will be shown\n")
    WORK_DIR <- getwd()
       
        observe({
         output$CDhistogram <- renderText({
           "CD histogram"
         })
          
          output$histPlotCD <- renderImage({
           req(input$whichHistCD)
           if (input$whichHistCD == "Normalized"){
             path <- paste0(WORK_DIR,"/7-metabolite_data_preprocessing/normalized/CD_histogram_norm.png")
             cat ("image path =",path,"\n")
           }
           if (input$whichHistCD == "Raw"){
             path <- paste0(WORK_DIR,"/7-metabolite_data_preprocessing/normalized/CD_histogram_raw.png")
             cat ("image path =",path,"\n")
           }
           
           list(src = path, contentType = 'image/png',width = "500px", height = "auto",
                alt = "This is alternate text")
           
         } ,deleteFile=FALSE)
      
          output$UChistogram <- renderText({
            "UC histogram"
          })
          
         output$histPlotUC <- renderImage({
           req(input$whichHistUC)
           if (input$whichHistUC == "Normalized"){
             path <- paste0(WORK_DIR,"/7-metabolite_data_preprocessing/normalized/UC_histogram_norm.png")
             cat ("image path =",path,"\n")
           }
           if (input$whichHistUC == "Raw"){
             path <- paste0(WORK_DIR,"/7-metabolite_data_preprocessing/normalized/UC_histogram_raw.png")
             cat ("image path =",path,"\n")
           }
           
           list(src = path, contentType = 'image/png',width = "500px", height = "auto",
                alt = "This is alternate text")
           
         } ,deleteFile=FALSE)

        })#observe 
    
  })#observeEvent
  
  
  
  
 
}#eof server