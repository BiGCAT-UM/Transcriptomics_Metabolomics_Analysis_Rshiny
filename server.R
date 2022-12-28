server = function(input, output,session) { 
  
  #to limit file size
  options(shiny.maxRequestSize = 500*1024^2)
  
  ################################################################
  
  # Transcriptomics
  
  # ################################################################
  # 
  # hideTab("tabs_trans", target = "filtering_trans")
  # hideTab("tabs_trans", target = "norm_trans")
  # hideTab("tabs_trans", target = "deg_trans")
  # hideTab("tabs_trans", target = "mapping_trans")
  # hideTab("tabs_trans", target = "pathway_trans")
  # hideTab("tabs_trans", target = "heatmap_trans")
  # hideTab("tabs_trans", target = "network_trans")
  # 
  # ################################################################
  # 
  # # Metabolomics
  # 
  # ################################################################
  # 
  hideTab("tabs_mets", target = "filtering_mets")
  #hideTab("tabs_mets", target = "norm_mets")
  #hideTab("tabs_mets", target = "stat_mets")
  #hideTab("tabs_mets", target = "pathway_mets")
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
          h3(strong("Histogram of Mean Expression Values (pre-filtering)"))
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
          h3(strong("Histogram of Mean Expression Values (post-filtering)"))
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
      text = "Data successfully filtered! Please wait for the figures to be rendered.",
      type = "success")
    
  })
  
  # Go the next step
  observeEvent(if (length(data_filtered())> 0){input$preprocess_NEXT}, {
    
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
      
      showModal(modalDialog(title = h4(strong("Normalization and Quality Control..."),
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
        text = "Normalization and QC successfully performed! 
        You can now view the generated QC plots.",
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
        
        
        showModal(modalDialog(title = h4(strong("Normalization and Quality Control..."),
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
          text = "Normalization and QC successfully performed! 
        You can now view the generated QC plots.",
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
        req(input$normalizedQC)
        req(input$colorQC)
        
        if (input$whichQCplot == "PCA"){
          if (input$normalizedQC == "Raw"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/PCAanalysis__biopsylocation2.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/PCAanalysis__disease2.png")
              cat ("image path =",path,"\n")
            }
          }
          if (input$normalizedQC == "Normalized"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/PCAanalysis__biopsylocation2.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/PCAanalysis__disease2.png")
              cat ("image path =",path,"\n")
            }
          }
        }
        if (input$whichQCplot == "Boxplot"){
          if (input$normalizedQC == "Raw"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/Boxplot__biopsylocation.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/Boxplot__disease.png")
              cat ("image path =",path,"\n")
            }
          }
          if (input$normalizedQC == "Normalized"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/Boxplot__biopsylocation.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/Boxplot__disease.png")
              cat ("image path =",path,"\n")
            }
          }
        }
        if (input$whichQCplot == "Heatmap"){
          if (input$normalizedQC == "Raw"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/Correlation__biopsylocation.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCraw/Correlation__disease.png")
              cat ("image path =",path,"\n")
            }
          }
          if (input$normalizedQC == "Normalized"){
            if (input$colorQC == "Location"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/Correlation__biopsylocation.png")
              cat ("image path =",path,"\n")
            }
            if (input$colorQC == "Disease"){
              path <- paste0(WORK_DIR,"/2-differential_gene_expression_analysis/QCnorm/Correlation__disease.png")
              cat ("image path =",path,"\n")
            }
          }
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
  FC_threshold <- eventReactive(input$DEGButton, {
    input$FCthreshold
  })
  
  # P value threshold
  P_threshold <- eventReactive(input$DEGButton, {
    input$pthreshold
  })
  
  # Read top table
  topTable_all <- eventReactive(input$DEGButton, {
    req(input$AdjOrRaw)
    topTable <- list()
    WORK_DIR <- getwd()
    
    # Loading message
    showModal(modalDialog(title = h4(strong("Statistical Analysis..."),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    # perform DE analysis
    DE_analysis(data_filtered()[[1]], data_filtered()[[2]], 0)
      
    # Read tables
    topTable[[1]] <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab"))
    topTable[[2]] <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab"))
    topTable[[3]]<- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab"))
    topTable[[4]] <- read.delim(paste0(WORK_DIR,"/2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab"))
    
    names(topTable) <- c("Ileum: CD vs non-IBD",
                         "Rectum: CD vs non-IBD",
                         "Ileum: UC vs non-IBD",
                         "Rectum: UC vs non-IBD")
    
    removeModal()
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "DEG analysis is finished!",
      type = "success")
    
    return(topTable)
  })
  
  # Filter top table
  topTable_filtered <- eventReactive(input$DEGButton, {
    req(input$AdjOrRaw)
    topTable <- list()
    WORK_DIR <- getwd()
    
    if (input$AdjOrRaw == "Raw") {
      temp <- topTable_all()[[1]]
      topTable[[1]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      temp <- topTable_all()[[2]]
      topTable[[2]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      temp <- topTable_all()[[3]]
      topTable[[3]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
      
      
      temp <- topTable_all()[[4]]      
      topTable[[4]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$pvalue < P_threshold()),]
    }
    if (input$AdjOrRaw == "Adjusted") {
      temp <- topTable_all()[[1]]
      topTable[[1]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      temp <- topTable_all()[[2]]
      topTable[[2]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      temp <- topTable_all()[[3]]
      topTable[[3]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
      
      temp <- topTable_all()[[4]]
      topTable[[4]] <- temp[(abs(temp$FoldChange) > FC_threshold()) &
                              (temp$padj < P_threshold()),]
    }

    names(topTable) <- c("Ileum: CD vs non-IBD",
                         "Rectum: CD vs non-IBD",
                         "Ileum: UC vs non-IBD",
                         "Rectum: UC vs non-IBD")
    
    return(topTable)
  })

  
  Comparison_DEG <- reactive({
    req(input$Comparison)
    return(input$Comparison)
  })
  
  RawOrAdj_DEG <- reactive({
    req(input$AdjOrRaw)
    return(input$AdjOrRaw)
  })
  
  # Return top table
  observe({
    output$topTable <- DT::renderDataTable({
      req(input$Comparison)
      output <- topTable_filtered()[[Comparison_DEG()]]
      output <- arrange(output, pvalue)
      output <- output[,c(1,2,3,4,7,8)]
      colnames(output) <- c("Gene Name", "Avg Expr", "log2 FC", "FC", "p-value", "adj. p-value")
      return(output)
    }, server=TRUE,
    options = list(pageLength = 5), rownames= FALSE)
    
  })

  output$VolcanoPlot <- renderPlot(NULL)
  
  observeEvent(input$DEGButton, {
    
    # Make plot
    output$VolcanoPlot <- renderPlot({
      req(input$AdjOrRaw)
      req(input$Comparison)
      topTable <- topTable_all()[[Comparison_DEG()]]
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
    
  })#eof observeEvent
  

  # Go the next step
  observeEvent(input$deg_NEXT, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "DEG analysis successfully completed! Now you can begin with the identifier mapping!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "mapping_trans")
    
    showTab("tabs_trans", target = "mapping_trans")
    
  })
  
  #***************************************************#
  # Identifier mapping
  #***************************************************#
  observeEvent(input$mappingButton,{
    # Loading message
    showModal(modalDialog(title = h4(strong("Identifier Mapping..."),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    mappingTranscriptomics(RawOrAdj_DEG())
    
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Identifiers have been successfully mapped!",
      type = "success")
    
    output$mappingTable <- DT::renderDataTable({
      output <- read.delim(paste0(work_DIR, "/3-identifier_mapping/IDMapping_",input$mappingDisease, ".tsv"))
      return(output)
    }, server=TRUE,
    options = list(pageLength = 10), rownames= FALSE)
    
  })
  
  # Go the next step
  observeEvent(input$mapping_NEXT, {
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "pathway_trans")
    
    showTab("tabs_trans", target = "pathway_trans")
    
  })
  #***************************************************#
  # Pathway analysis
  #***************************************************#

  observeEvent(input$pathwayButton,{
    
    showModal(modalDialog(title = h4(strong("Pathway Analysis..."),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    # Perform pathway analysis
    pathwayAnalysisTranscriptomics(log2(FC_threshold()), P_threshold(), 0.05,0.02)
    
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Pathway analysis has been performed successfully!",
      type = "success")
    
    # pathway table
    output$pathwayTable <- DT::renderDataTable({
      req(input$pathwayComparison)
      if (input$pathwayComparison == "Ileum: CD vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_CD_ileum.tsv"))
      }
      if (input$pathwayComparison == "Rectum: CD vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_CD_rectum.tsv"))
      }
      if (input$pathwayComparison == "Ileum: UC vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_UC_ileum.tsv"))
      }
      if (input$pathwayComparison == "Rectum: UC vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_UC_rectum.tsv"))
      }
      
      output <- output[,1:7]
      colnames(output) <- c("ID", "Description", "Gene ratio", "Bg ratio", "p-value","adj. p-value",  "q-value")
      return(output)
    }, server=TRUE,
    options = list(pageLength = 5), rownames= FALSE)
    
    
    output$pathwayPlot <- renderPlot({
      req(input$pathwayComparison)
      if (input$pathwayComparison == "Ileum: CD vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_CD_ileum.tsv"))
      }
      if (input$pathwayComparison == "Rectum: CD vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_CD_rectum.tsv"))
      }
      if (input$pathwayComparison == "Ileum: UC vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_UC_ileum.tsv"))
      }
      if (input$pathwayComparison == "Rectum: UC vs non-IBD"){
        output <- read.delim(paste0(work_DIR,"/4-pathway_analysis/enrichResults_ORA_UC_rectum.tsv"))
      }
      
      output <- arrange(output, pvalue)
      p <- ggplot(output[1:5,]) +
        geom_segment(aes(y = 0, yend = -log10(pvalue), x = reorder(Description,-1*pvalue), xend = Description)) +
        geom_point(aes(y = -log10(pvalue), x = reorder(Description, -1*pvalue), 
                       size = Count, color = -log10(pvalue))) +
        coord_flip() +
        labs(size = "Set Size", color = "-log10 p-value", x = "", y = "-log10 p-value") +
        scale_color_viridis_c() +
        theme_minimal() +
        theme(axis.text = element_text(size = 20))
      
      return(p)
    })
    

  })

  # Go the next step
  observeEvent(input$pathway_NEXT, {
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "heatmap_trans")
    
    showTab("tabs_trans", target = "heatmap_trans")
    
  })
  
  #***************************************************#
  # Create heatmap
  #***************************************************#
  
  p_threshold_pathway <- eventReactive(input$heatmapButton, {
    input$p_threshold_pathway
  })
  
  q_threshold_pathway <- eventReactive(input$heatmapButton, {
    input$q_threshold_pathway
  })
  
  observeEvent(input$heatmapButton, {
    
    showModal(modalDialog(title = h4(strong("Heatmap..."),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    createHeatmap(p_threshold_pathway(),q_threshold_pathway())
    
    removeModal()
    
    output$HeatmapPlot<- renderImage({
      path <- paste0(work_DIR,"/5-create_heatmap/heatmap_log10_large.png")
      list(src = path, contentType = 'image/png',width = "800px", height = "auto",
           alt = "This is alternate text")
      
    }, deleteFile=FALSE)
    
  })
  
  # Go the next step
  observeEvent(input$heatmap_NEXT, {
    
    updateTabsetPanel(session, "tabs_trans",
                      selected = "network_trans")
    
    showTab("tabs_trans", target = "network_trans")
    
  })
  
  #***************************************************#
  # Network analysis
  #***************************************************#
  
  observeEvent(input$networkButton, {
    
    showModal(modalDialog(title = h4(strong("Network Analysis"),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    networkAnalysis()
    
    removeModal()
    
    output$NetworkPlot<- renderImage({
      wp.hs.gmt <- "wikipathways-20220510-gmt-Homo_sapiens.gmt"
      path <- paste0("6-network_analysis/PPI_Pathway_Network_", wp.hs.gmt, input$location_network,".png")
      list(src = path, contentType = 'image/png',width = "800px", height = "auto",
           alt = "This is alternate text")
      
    }, deleteFile=FALSE)
    
  })
  
  
  
  
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
  
  shapiroResults <- eventReactive (input$normButton,{ 
    
    #call function     
    shapiroResults <- normalizeMets(selectedMethod())
  
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Normalization and QC plots done. 
        You can find all QC plots in 7-metabolite_data_preprocessing",
      type = "success")
  return (shapiroResults)  
    
  })#eventReactive
  
  observeEvent (length(shapiroResults())>0,{
     
      output$CDhistogram <- renderUI({
        tagList(
          h3(strong("CD histogram"))
        )
      })
      
      output$CD_shapiro_result <- renderText({
        if(shapiroResults()[[1]] == TRUE){
          paste0("Data after transformation seems to follow a normal distribution for more then 80% of your data")}
        else{
          "Advised to select a different data transformation procedure"
        }
      })#renderText
   
      WORK_DIR <- getwd()
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
    
      output$CDhistogram <- renderUI({
        tagList(
          h3(strong("CD histogram"))
        )
      })
      
      output$UC_shapiro_result <- renderText({
        if(shapiroResults()[[2]] == TRUE){
          paste0("Data after transformation seems to follow a normal distribution for more then 80% of your data")}
        else{
          "Advised to select a different data transformation procedure"
        }
      })#renderText
      
      WORK_DIR <- getwd()
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
      
    
  })#observeEvent
  
  # Go the next step
  observeEvent(input$metNorm_NEXT, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Normalization successfully completed! Now you can continue with statistical analysis!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_mets",
                      selected = "stat_mets")
    
    showTab("tabs_mets", target = "stat_mets")
    
  })
 
  
   
  #***************************************************#
  # Statistical Analysis
  #***************************************************#
  
  FC_met <- eventReactive(input$statButton,{
    input$FCthresholdMet
  })
  pvalue_met <- eventReactive(input$statButton,{
    input$pthresholdMet
  })
  
  
  #to get all results from statistical analysis
  resTable <- eventReactive(input$statButton,{
    
    allResults <- list()
   
    # Loading message
    showModal(modalDialog(title = h4(strong("Metabolomics statistical analysis started"),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    
    mSet_transformedCD <- read.csv("7-metabolite_data_preprocessing/normalized/CD_norm_data.csv", na.strings=c("", "NA"))
    mSet_transformedUC <- read.csv("7-metabolite_data_preprocessing/normalized/UC_norm_data.csv", na.strings=c("", "NA"))
    
    statAnalysisMets(mSet_transformedCD, "CD", selectedMethod(), FC_met(),pvalue_met())
    statAnalysisMets(mSet_transformedUC, "UC", selectedMethod(), FC_met(),pvalue_met())
 
    allResults [[1]] <- read.csv("8-significantly_changed_metabolites_analysis/mbxData_CD.csv")
    allResults [[2]] <- read.csv("8-significantly_changed_metabolites_analysis/mbxData_UC.csv")
    
    names(allResults) <- c("CD vs non-IBD",
                        "UC vs non-IBD")
    
    removeModal()
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Statistical analysis finished!",
      type = "success")
    
    
    return (allResults)
  })
  
  filteredResults <- eventReactive(input$statButton, {
    
    filtered <- list()
    #WORK_DIR <- getwd()
    
    #if (input$AdjOrRaw == "Raw") {
      temp <- resTable()[[1]]
      filtered[[1]] <- temp[(abs(temp$foldchange_disorder) > FC_met()) &
                              (temp$p_values_disorder < pvalue_met()),]
#browser()
      temp <- resTable()[[2]]
      filtered[[2]] <- temp[(abs(temp$foldchange_disorder) > FC_met()) &
                              (temp$p_values_disorder < pvalue_met()),]
      
    #}
      names(filtered) <- c("CD vs non-IBD",
                          "UC vs non-IBD")
      return (filtered)
      
  })
  
  compPairMet <- reactive({
    req(input$metCompPair)
    return(input$metCompPair)
  })
  
  # Show filtered result table for selected comparison
  observe({
    
    output$metResTable <- DT::renderDataTable({
      req(input$metCompPair)
      output <- filteredResults()[[compPairMet()]]
      # output <- arrange(output, pvalue)
      # output <- output[,c(1,2,3,4,7,8)]
      # colnames(output) <- c("Gene Name", "Avg Expr", "log2 FC", "FC", "p-value", "adj. p-value")
      return(output)
    }, server=TRUE,
    options = list(pageLength = 5), rownames= FALSE)
    
  })
  
  output$metVolcanoPlot <- renderPlot(NULL)
  
  #show volcano plot for selected comparison
  
  
  observeEvent(input$statButton,{
    
      output$metVolcanoPlot <- renderImage({
        WORK_DIR <- getwd()
        req(input$metCompPair)
        if (compPairMet() == "CD vs non-IBD"){
         # path <- paste0(WORK_DIR,"/8-significantly_changed_metabolites_analysis/CD_relevant_labels_VolcanoPlot_absLogFC_0.58_pValue_0.05.png")
          path <- paste0(WORK_DIR,"/8-significantly_changed_metabolites_analysis/CD_relevant_labels_VolcanoPlot.png")
          cat ("image path =",path,"\n")
        }
        if (compPairMet() == "UC vs non-IBD"){
         # path <- paste0(WORK_DIR,"/8-significantly_changed_metabolites_analysis/UC_relevant_labels_VolcanoPlot_absLogFC_0.58_pValue_0.05.png")
          path <- paste0(WORK_DIR,"/8-significantly_changed_metabolites_analysis/UC_relevant_labels_VolcanoPlot.png")
          cat ("image path =",path,"\n")
        }
        
        list(src = path, contentType = 'image/png',width = "500px", height = "auto",
             alt = "This is alternate text")
        
      } ,deleteFile=FALSE)
  })
  
  # Go the next step
  observeEvent(input$stat_NEXT, {
    
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Statistical analysis successfully completed! 
      Now you can continue with pathway analysis!",
      type = "success")
    
    updateTabsetPanel(session, "tabs_mets",
                      selected = "pathway_mets")
    
    showTab("tabs_mets", target = "pathway_mets")
    
  })
  
  #***************************************************#
  # Pathway Analysis
  #***************************************************#
  
  observeEvent(input$pathwayButtonMet,{
    
    
    
    
  })#observeEvent
  
 
}#eof server