work_DIR <- getwd()
# check if BioCmanager libraries are already installed > otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"baySeq" %in% installed.packages()) BiocManager::install("baySeq")
if(!"DESeq2" %in% installed.packages()) BiocManager::install("DESeq2")
if(!"edgeR" %in% installed.packages()) BiocManager::install("edgeR")
if(!"bioDist" %in% installed.packages()) BiocManager::install("bioDist")
if(!"biomaRt" %in% installed.packages()) BiocManager::install("biomaRt")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"magrittr" %in% installed.packages()) BiocManager::install("magrittr")
if(!"EnhancedVolcano" %in% installed.packages()) BiocManager::install("EnhancedVolcano")
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2")}
if(!"limma" %in% installed.packages()){install.packages("limma")}
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if(!"readxl" %in% installed.packages()) BiocManager::install("readxl")
if(!"downloader" %in% installed.packages()) BiocManager::install("downloader")
if(!"R.utils" %in% installed.packages()){install.packages("R.utils")}
if(!"pheatmap" %in% installed.packages()) BiocManager::install("pheatmap")
if(!"BridgeDbR" %in% installed.packages()) BiocManager::install("BridgeDbR")

library(rstudioapi)
library(readxl)
library(dplyr)
library(magrittr)
library(edgeR)#avgLogCPM function
library(ggplot2)
library(baySeq)
library(DESeq2)
library(bioDist)
library(biomaRt)
library(EnhancedVolcano)
library(limma)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(shinyBS)
library(shinycssloaders)
library(downloader)
library(R.utils)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(RCy3)
library(BridgeDbR)
library(rJava)

source("functions_ArrayAnalysis_v2.R")

#Normalization and QC plot generation
normalize_QCplots <- function (sampleLabels, htxCount){
  
  sampleLabels$disease <- relevel(factor(sampleLabels$disease),ref="nonIBD")
  #add an experimental group variable to sampleLabels
  sampleLabels$group <- as.factor(paste(sampleLabels$disease,sampleLabels$biopsy_location,sep="_"))
  
  # First create a DESeqDataSet object
  #(non-intercept) statistical model based on the disease and biopsy_location, group column represent both of them
  dds <- DESeqDataSetFromMatrix(countData = htxCount, colData=sampleLabels, design= ~0 + group)
  
  #estimate the size factors
  #To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us.
  dds <- estimateSizeFactors(dds)
  #normalize the data (here for Quality Control(QC) plotting)
  #QC plotting is optional
  norm <- counts(dds,normalize=TRUE)
  #create a 2logged data for original object (here for QC plotting)
  datlog <- log(htxCount+1,2)
  #create a 2logged norm object (here for QC plotting)
  normlog <- log(norm+1,2)
  #for QC remove genes that have not been measured in any sample in the experiment
  datlogQC <- datlog[rowSums(datlog)!=0,]
  normlogQC <- normlog[rowSums(normlog)!=0,]
  
  cat ("Normalization done\n ")
  
  WORK.DIR <- getwd()  
  setwd(paste0(WORK.DIR,"/2-differential_gene_expression_analysis")) 
  #write.table(normlogQC, "normExpression.csv", sep=",",quote=FALSE,row.names = TRUE)

  #create QC plots for raw data, colored by different variables
  factors <- c("disease","biopsy_location","group")
  if(!dir.exists("QCraw")){dir.create("QCraw")}
  setwd("QCraw")
  png("sizefactors.png")
  plot(sizeFactors(dds),type='h',lwd=5,ylim=c(0,max(sizeFactors(dds))),col="darkblue")
  dev.off()
  createQCPlots(datlogQC, factors, Table=sampleLabels, normMeth="", postfix="")
  setwd("..")#go back to main directory
  
  cat ("QC plot for raw data done\n ")
  
  #create QC plots for normalized data colored by different variables
  if(!dir.exists("QCnorm")){dir.create("QCnorm")}
  setwd("QCnorm")
  createQCPlots(normlogQC, factors, Table=sampleLabels, normMeth="DESeq", postfix="")
  setwd("..")#go back to main dir
  
  cat ("QC plot for normalized data done\n ")
  
  #go back main directory 
  setwd('..')
}

#Removing outlier samples
removeOutliers <- function (sampleLabels, htxCount,outliers){

  #remove selected outliers from both orj data and sample Labels
  htxCount <- htxCount[,-match(outliers,colnames(htxCount))]
  sampleLabels <- sampleLabels[-match(outliers,rownames(sampleLabels)),]
  
  
  return(list((sampleLabels),(htxCount)))
}

#Filtering sample and genes with all zero values
sample_gene_filtering<-function(htxMeta,htxOrj){
  
  # #to get current working directory
  # WORK.DIR <- getwd()  
  # 
  # # set working environment to the location where current source file is saved into.
  # setwd(paste0(WORK.DIR,"/1-data_preprocessing")) 
  
  #filter out samples by data type as host_transcriptomics
  htxMeta <- htxMeta  %>% filter(htxMeta$data_type == "host_transcriptomics")
  
  #filter out data by biopsy location, include CD, UC and nonIBD samples from ileum and rectum location 
  htxMeta <-htxMeta  %>% filter (  (htxMeta$diagnosis == "CD" & htxMeta$biopsy_location=="Ileum") 
                                   | (htxMeta$diagnosis == "CD" & htxMeta$biopsy_location=="Rectum")
                                   | (htxMeta$diagnosis == "UC" & htxMeta$biopsy_location=="Ileum") 
                                   | (htxMeta$diagnosis == "UC" & htxMeta$biopsy_location=="Rectum") 
                                   | (htxMeta$diagnosis == "nonIBD" & htxMeta$biopsy_location=="Rectum") 
                                   | (htxMeta$diagnosis == "nonIBD" & htxMeta$biopsy_location=="Ileum") 
  )
  #filter out samples by visit_num=1
  htxMeta <-htxMeta  %>% filter(htxMeta$visit_num == "1")
  #filter out unused columns
  htxMeta <- htxMeta %>% dplyr::select(External.ID,Participant.ID,biopsy_location,diagnosis)
  #Order htxMeta data based on external ID to match samples with htx count correctly
  htxMeta<- htxMeta[order(htxMeta$External.ID),]#order htxMeta by external ID
  
  #Convert sample names to upper (some of them are in lower case)
  colnames(htxOrj)<-toupper(colnames(htxOrj))
  #htx count data is filtered based on column names in htxMeta
  names.use <- names(htxOrj)[(names(htxOrj) %in% htxMeta$External.ID)]
  #filter out htxOrj based on names.use and create a new htxCount
  htxCount <- htxOrj[, names.use]
  #htxCount data are ordered based on column names to match samples between htxCount and sampleLabels
  htxCount <- htxCount[,order(names(htxCount))]
  
  #checking which samples have all zero values across all genes
  #these sample should be removed otherwise there will be a problem when calculating estimate size factors
  idx <- which(colSums(htxCount) == 0)
  #CSMDRVXI MSM719ME  are samples which has all zero values for all genes, so we remove them
  htxCount <- htxCount[ , -idx]
  
  #removing same samples from sample labels metadata
  htxMeta <- htxMeta[-idx , ]
  #Set column one as rownames
  rownames(htxMeta) <- htxMeta[,1]
  htxMeta <- htxMeta[,c(-1,-2)]
  #add column names 
  colnames(htxMeta) <- c( "biopsy_location","disease")
  #check whether sample names are in same order
  #all(colnames(htxCount) == rownames(sampleLabels2))
  
  #remove genes which have all zero values across all samples then start DE analysis
  nonzero <- rowSums(htxCount) > 0
  htxCount %<>% .[nonzero,]
  
# setwd('..')

return(list((htxMeta),(htxCount)))
       
}

#Filtering applying CPM method
cpm_filter<-function(htxMeta, htxCount, filter_threshold){
#selected threshold should be converted to numeric value
filter_threshold <- as.numeric(filter_threshold)
#aveLogCPM function computes average log2 counts-per-million for each row of counts.
#the below function is similar to log2(rowMeans(cpm(y, ...)))  
mean_log_cpm = aveLogCPM(htxCount)

cat("selected threshold:", filter_threshold,"\n")

# We plot the distribution of average log2 CPM values to verify that our chosen presence threshold is appropriate. 
#The distribution is expected to be bi modal, with a low-abundance peak representing non-expressed genes and a high-abundance peak representing expressed genes. The chosen threshold should separate the two peaks of the bi modal distribution. 

#jpeg(file="avgLogCpmDist.jpeg")#if you want to save the histogram uncomment the following command  
histogram_tmp <- ggplot() + aes(x=mean_log_cpm) +
  geom_histogram(binwidth=0.2, fill = "grey", color = "black") +
  geom_vline(xintercept = filter_threshold, color = "red", linewidth = 1.5, linetype = "dashed") +
  theme_classic() +
  xlab("Mean log CPM") +
  ylab("Count")
  #ggtitle("Histogram of mean expression values")
#dev.off()#to save the plot to the file
#Having chosen our threshold, lets pick the subset of genes whose average expression passes that threshold.
keep_genes <- mean_log_cpm >= filter_threshold 
htxCount <- htxCount[keep_genes,]
cat("\nortada\n")
#sample distribution based on biopsy locations
ileum = nrow(htxMeta[htxMeta$biopsy_location=="Ileum",])
rectum = nrow(htxMeta[htxMeta$biopsy_location=="Rectum",])
cat ("Number of samples in ileum:", ileum ,"\nNumber of samples in rectum:",rectum)

#create output folder if it doesn't exist
#if(!dir.exists("output")) dir.create("output")

#Write all the generated data into the related output files 
write.table(htxCount, "1-data_preprocessing/htxCount.csv", sep=",",quote=FALSE, row.names = TRUE )
write.table(htxMeta, "1-data_preprocessing/sampleLabels.csv", sep=",",quote=FALSE,row.names = TRUE)

cat("\nPreprocessing is finished, results are saved to output folder.")

return (histogram_tmp)

}

#Filtering applying CPM method
cpm_filter_output <- function(htxMeta, htxCount, filter_threshold){
  
  #selected threshold should be converted to numeric value
  filter_threshold <- as.numeric(filter_threshold)
  
  #aveLogCPM function computes average log2 counts-per-million for each row of counts.
  #the below function is similar to log2(rowMeans(cpm(y, ...)))  
  mean_log_cpm <- aveLogCPM(htxCount)
  
  #Having chosen our threshold, lets pick the subset of genes whose average expression passes that threshold.
  keep_genes <- mean_log_cpm >= filter_threshold 
  htxCount <- htxCount[keep_genes,]
  output <- list(htxMeta, htxCount)
  return(output)
}

#Differential Gene Expression Analysis   
DE_analysis <- function (sampleLabels, htxCount, FC_threshold){

WORK.DIR <- getwd()
#browser()
#this line of code can be removed
#htxCount <- read.csv(paste0(WORK.DIR,"/1-data_preprocessing/htxCount.csv"))
#sampleLabels <- read.csv(paste0(WORK.DIR,"/1-data_preprocessing/sampleLabels.csv"),row.names = 1 )

setwd(paste0(WORK.DIR,"/2-differential_gene_expression_analysis"))
sampleLabels$disease <- relevel(factor(sampleLabels$disease),ref="nonIBD")
#add an experimental group variable to sampleLabels
sampleLabels$group <- as.factor(paste(sampleLabels$disease,sampleLabels$biopsy_location,sep="_"))

dds <- DESeqDataSetFromMatrix(countData = htxCount, colData=DataFrame(sampleLabels), design= ~0 + group)
dds <- estimateSizeFactors(dds)

################################################################################
#######################statistical modelling####################################
################################################################################
#set directory for stat output
if(!dir.exists("statsmodel")) dir.create("statsmodel")
setwd("statsmodel")
#run differential analysis
dds <- DESeq(dds)
cont.matrix <- makeContrasts(
  #CD disease on ileum and rectum
  CD_Ileum_vs_nonIBD_Ileum   = groupCD_Ileum - groupnonIBD_Ileum,
  CD_Rectum_vs_nonIBD_Rectum = groupCD_Rectum - groupnonIBD_Rectum,
  #UC disease on ileum and rectum
  UC_Ileum_vs_nonIBD_Ileum   = groupUC_Ileum - groupnonIBD_Ileum,
  UC_Rectum_vs_nonIBD_Rectum = groupUC_Rectum - groupnonIBD_Rectum,
  #UC and CD disease comparison
  UC_Ileum_vs_CD_Illeum     = groupUC_Ileum - groupCD_Ileum,
  UC_Rectum_vs_CD_Rectum    = groupUC_Rectum - groupCD_Rectum,
  #biopsy location comparisons
  CD_Rectum_vs_CD_Ileum     = (groupCD_Rectum - groupnonIBD_Rectum) - (groupCD_Ileum - groupnonIBD_Ileum),
  UC_Rectum_vs_UC_Ileum     = (groupUC_Rectum - groupnonIBD_Rectum) - (groupUC_Ileum - groupnonIBD_Ileum),
  levels = resultsNames(dds)
)
#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
#the function results() is called from within the saveStatOutputDESeq2 function to compute the contrasts
 files <- saveStatOutputDESeq2(cont.matrix,dds,postfix="",annotation=NULL)
# #create summary table of the contrast results
#
# cat ("FC degeri", FC_threshold,"\n")
createPvalTab(files,postfix="",namePVal="pvalue",nameAdjPVal="padj",nameFC="FoldChange",nameLogFC="log2FoldChange",html=FALSE, FC_threshold)

# 
#   WORK.DIR <- getwd()
#  # setwd(paste0(WORK.DIR,"/2-differential_gene_expression_analysis"))
#   readPath <- paste0(WORK.DIR,"/2-differential_gene_expression_analysis/statsmodel/Summary_tables.tab")
#   summaryTab <- readLines(readPath)
#   
#   for (i in 2:length(summaryTab)){
#     strsplit(summaryTab[i], split = "\t")
#     
#   }
  
}

#
showFileList <- function(){
  
  WORK.DIR <- getwd()
  #path of read folder  
  readFilePath <- paste(WORK.DIR,"2-differential_gene_expression_analysis/statsmodel",sep="/")
  #read all files with .tab extension under statsmodel folder
  allFiles <- list.files(path=readFilePath, pattern = ".tab", all.files=TRUE,full.names=TRUE)
  
  fileList = list()
  for(i in 1:length(allFiles)) {
    splitted <- strsplit(allFiles[i], split = "/")[[1]]
    #if last value starts with "table" then we will keep it
    if(startsWith(splitted[[length(splitted)]],"table")){
      cat(splitted[[length(splitted)]],"\n")
      fileList <- append(fileList, splitted[[length(splitted)]][[1]])
    }
  }
  #browser()
  # convert file list to data frame to be shown in UI
  fileList <- data.frame(matrix(unlist(fileList)),stringsAsFactors=FALSE)
  
  #fileNames and full paths are send to back
  return (fileList)
}

#to show volcano plot of selected comparison pair
volcanoPlots <- function (selected = ""){
  
 #browser()
 cat (selected)
  # ## Enhanced Volcano Plots Visualization
  # #Differential gene expression analysis results is visualized by volcano plots
  # WORK.DIR <- getwd()
  # readFilePath <- paste(WORK.DIR,"statsmodel",sep="/")
  # # create an empty list
  # plot_list = list()
  # 
  #  for(i in 1:length(files))
  #  {
  #     #read file
  #     splitted <- strsplit(files[i], split = "_")[[1]]
  #     title <- splitted [2:6]
  #     title <- paste (splitted[2],splitted[3],splitted[4],splitted[5],splitted[6],sep=" ")
  #     title <- (strsplit(title, split = "\\.")[[1]])[1]
  #     tab <- read.delim(paste(readFilePath, files[i],sep="/"),header=TRUE,as.is=TRUE)
  #     
  #     p <- EnhancedVolcano(tab , lab = tab$X, labSize = 3, title = title,
  #                          x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, FCcutoff = 0.58)
  #     plot_list[[i]] = p
  #  }
  # 
  # #path of the output folder
  # outFolder <- paste(WORK.DIR,"volcano_plots",sep="/")
  # #create folder if doesn't exist
  # if(!dir.exists(outFolder)) dir.create(outFolder)
  # 
  # #for(i in 1:length(files)) {
  # file_name = paste(outFolder,"/",files[i],".png",sep="")
  # png(file_name)
  # print(plot_list[[i]])
  # dev.off()
  # # }
  # #go back to main folder
  # setwd('..')

#  return (files)
}
#identifier mapping of transcriptomics data 
###################Identifier Mapping for DE analysed data#############################
mappingTranscriptomics <- function (RawOrAdj){
  
  #input data is all files starts with "table" under the folder 2-differential_gene_expression_analysis/statsmodel/
  #output files will be in "3-identifier_mapping" folder with the same structure IDMapping_CD OR IDMapping_UC
  
  # Read top tables
  setwd(work_DIR)
  dataset1 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab")
  dataset2 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab")
  dataset3 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab")
  dataset4 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab")
  
  disease <- c("CD", "UC")
  
  for (d in disease){
    if (RawOrAdj == "Raw"){
      if (d == "CD") {
        #filter out  unused columns, we select geneSymbol, log2FC and pvalue
        dataset_ileum<- subset( dataset3, select = c(1,3,7))
        dataset_rectum<- subset( dataset4, select = c(1,3,7))
        print("Selected disorder is Crohn's disease")
      }else if(d == "UC"){ 
        #filter out  unused columns, we select geneSymbol, log2FC and pvalue
        dataset_ileum<- subset( dataset1, select = c(1,3,7))
        dataset_rectum<- subset( dataset2, select = c(1,3,7))
        print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
        }
    }
    if (RawOrAdj == "Adjusted"){
      if (d == "CD") {
        #filter out  unused columns, we select geneSymbol, log2FC and pvalue
        dataset_ileum<- subset( dataset3, select = c(1,3,7))
        dataset_rectum<- subset( dataset4, select = c(1,3,7))
        print("Selected disorder is Crohn's disease")
      }else if(d == "UC"){ 
        #filter out  unused columns, we select geneSymbol, log2FC and pvalue
        dataset_ileum<- subset( dataset1, select = c(1,3,7))
        dataset_rectum<- subset( dataset2, select = c(1,3,7))
        print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
        }
    }
    
    #merge two dataset of two locations into one data 
    dataset <- merge(dataset_ileum, dataset_rectum,by.x="X", by.y="X",sort = TRUE, all.x = TRUE, all.y = TRUE)
    
    #change column names
    colnames(dataset) <- c("GeneSymbol","log2FC_ileum","pvalue_ileum","log2FC_rectum","pvalue_rectum")
    
    ####################
    # HGNC to ENTREZ ID
    ####################
    hs <- org.Hs.eg.db #This object is a simple mapping of Entrez Gene identifier
    entrezID <- AnnotationDbi::select(hs, keys = dataset$GeneSymbol,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "SYMBOL")
    #filter out double gene symbols
    entrezID <- entrezID %>% distinct (entrezID$SYMBOL, .keep_all = TRUE)
    # add entrezIDs for each gene symbol in the dataset
    dataset <- cbind(entrezID$ENTREZID,dataset)
    #change column name
    colnames(dataset)[1] = "ENTREZ.ID"
    
    ####################
    # HGNC to ENSEMBL ID
    ####################
    hs <- org.Hs.eg.db #This object is a simple mapping of Entrez Gene identifier
    
    ensemblID <- AnnotationDbi::select(hs, keys = dataset$GeneSymbol,
                                       columns = c("ENSEMBL", "SYMBOL"),
                                       keytype = "SYMBOL")
    #filter out double gene symbols
    ensemblID <- ensemblID %>% distinct (ensemblID$SYMBOL, .keep_all = TRUE)
    
    # add entrezIDs for each gene symbol in the dataset
    dataset <- cbind(ensemblID$ENSEMBL,dataset)
    #change column name
    
    colnames(dataset)[1] = "Ensembl.ID"
    
    # Write output
    if(!dir.exists("3-identifier_mapping")){dir.create("3-identifier_mapping")}
    write.table(dataset, file=paste0("3-identifier_mapping/IDMapping_",d, ".tsv"),
                sep = "\t" ,quote = FALSE, row.names = FALSE)
  }
  
}

# Pathway analysis
pathwayAnalysisTranscriptomics <- function(logFCtheshold, Pthreshold, Pthreshold_pathway, Qthreshold_pathway){
  setwd(work_DIR)
  for (disorder in c("CD", "UC")){
    # Load data
    if (disorder == "CD") {
      dataset_CD <- read.delim("3-identifier_mapping/IDMapping_CD.tsv")
      #filter out  unused columns, we select Entrez.ID, log2FC and pvalue, remove NA values: #background genes to be used in enrichment analysis
      dataset <- na.omit(subset( dataset_CD, select = c(3,2,4:7)))
      print("Selected disorder is Crohn's disease")
    }else if(disorder == "UC"){ 
      dataset_UC <- read.delim("3-identifier_mapping/IDMapping_UC.tsv")
      #filter out  unused columns, we select Entrez.ID, log2FC and pvalue
      dataset <- na.omit(subset( dataset_UC, select = c(3,2,4:7)))
      print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
      }
    
    if(!dir.exists("4-pathway_analysis")) dir.create("4-pathway_analysis")
    #we will use selection criteria as Fold change=1.5,log2FC=0.58 and p.value < 0.05
    #for ileum location
    up.genes.ileum   <- dataset[dataset$log2FC_ileum >= logFCtheshold & dataset$pvalue_ileum < Pthreshold, 2] 
    down.genes.ileum <- dataset[dataset$log2FC_ileum <= -1*logFCtheshold & dataset$pvalue_ileum < Pthreshold, 2] 
    deg.ileum <- unique(dataset[(!is.na(dataset$ENTREZ.ID)) & (!is.na(dataset$pvalue_ileum)) & (dataset$pvalue_ileum < Pthreshold) & (abs(dataset$log2FC_ileum) > logFCtheshold),c(1:4)])
    write.table(deg.ileum,file = paste0("4-pathway_analysis/DEGs_",disorder,"_ileum.tsv"),sep="\t", quote=FALSE, row.names = FALSE)
    #for rectum location
    up.genes.rectum   <- dataset[dataset$log2FC_rectum >= logFCtheshold & dataset$pvalue_rectum < Pthreshold, 2] 
    down.genes.rectum <- dataset[dataset$log2FC_rectum <= -1*logFCtheshold & dataset$pvalue_rectum < Pthreshold, 2] 
    deg.rectum <- unique(dataset[!is.na(dataset$ENTREZ.ID) & !is.na(dataset$pvalue_rectum) & dataset$pvalue_rectum < Pthreshold & abs(dataset$log2FC_rectum) > logFCtheshold,c(1,2,5,6)])
    write.table(deg.rectum, file=paste0("4-pathway_analysis/DEGs_",disorder,"_rectum.tsv"),sep="\t", quote=FALSE, row.names = FALSE)
    
    pathway_data <- "local" #Options: local, new
    if (pathway_data == "local") {
      wp.hs.gmt <-list.files(work_DIR, pattern="wikipathways", full.names=FALSE)
      paste0("Using local file, from: ", wp.hs.gmt )
    }else if(pathway_data == "new"){ 
      #below code should be performed first to handle the ssl certificate error while downloading pathways 
      options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
      #downloading latest pathway gmt files for human 
      wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
      paste0("Using new data, from: ", wp.hs.gmt)}else{print("Pathway data type not recognized")
      }
    
    wp2gene   <- rWikiPathways::readPathwayGMT(wp.hs.gmt)
    wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
    ewp.ileum <- clusterProfiler::enricher(
      deg.ileum$ENTREZ.ID,# a vector of gene IDs
      universe = as.character(dataset$ENTREZ.ID), #background genes to be used in enrichment analysis
      pAdjustMethod = "fdr",#you can change this to hochberg, bonferronni or none etc.
      pvalueCutoff = 1, #adjusted pvalue cutoff on enrichment tests to report, we set it a wider criteria then we will filter out
      #results based on padjust and qvalue in next section which is enrichment result visualization
      qvalueCutoff = 1, #qvalue cutoff on enrichment tests to report as significant, 
      #multiple hypothesis testing
      TERM2GENE = wpid2gene, #user input annotation of TERM TO GENE mapping
      TERM2NAME = wpid2name) #user input of TERM TO NAME mapping
    ewp.ileum.res <- as.data.frame(ewp.ileum) 
    #Count all significant pathways that have p.adjust value lower than 0.05 and qvalue<0.02
    ileum.sign <- ewp.ileum.res[(ewp.ileum.res$p.adjust<Pthreshold_pathway)&(ewp.ileum.res$qvalue<Qthreshold_pathway),]
    #interpretation of the output: 
    #     BgRatio   = (number of genes measured in the current pathway) / (number of genes measured in all pathways)
    #     geneRatio = (number of DEGs in the current pathway) / (total number of DEGs in all pathways)
    ##Print location:
    paste0("Pathways enrichment results for disorder: ", disorder , ", location: ILEUM")
    # number of genes measured in all pathways
    paste0("The number of genes measured in all pathways is: ", length(ewp.ileum@universe))
    # number of DEGs in all pathways
    paste0("The number of DEGs measured in all pathways is: ", length(deg.ileum$ENTREZ.ID[deg.ileum$ENTREZ.ID %in% unique(wp2gene$gene)]))
    #number of enriched pathways
    paste0("The number of enriched pathways is: ", num.pathways.ileum <- dim(ewp.ileum.res)[1])
    #number of significantly enriched pathways
    paste0("The number of significantly enriched pathways is: ", num.pathways.ileum.sign <- dim(ileum.sign)[1])
    #exporting results to the file
    write.table(ewp.ileum.res, file=paste0("4-pathway_analysis/enrichResults_ORA_",disorder,"_ileum.tsv"),
                sep = "\t" ,quote = FALSE, row.names = FALSE)
    ##################RECTUM location#######################
    ewp.rectum <- clusterProfiler::enricher(
      deg.rectum$ENTREZ.ID, # a vector of gene IDs; options
      universe = as.character(dataset$ENTREZ.ID), #background genes to be used in enrichment analysis
      pAdjustMethod = "fdr",#you can change it as BH, bonferronni etc.
      pvalueCutoff = 1, #padjust cutoff
      qvalueCutoff = 1, #q value cutoff 
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name)
    ewp.rectum.res <- as.data.frame(ewp.rectum) 
    #Count all significant pathways that have p.adjust value lower than 0.05 and qvalue<0.02
    rectum.sign <- ewp.rectum.res[(ewp.rectum.res$p.adjust<Pthreshold_pathway)&(ewp.rectum.res$qvalue<Qthreshold_pathway),]
    ##Print location:
    paste0("Pathways enrichment results for disorder: ", disorder , ", location: RECTUM")
    # number of genes measured in all pathways
    paste0("The number of genes measured in all pathways is: ", length(ewp.rectum@universe))
    # number of DEGs in all pathways
    paste0("The number of DEGs measured in all pathways is: ", length(deg.rectum$ENTREZ.ID[deg.rectum$ENTREZ.ID %in% unique(wp2gene$gene)]))
    #number of enriched pathways
    paste0("The number of enriched pathways is: ", num.pathways.rectum <- dim(ewp.rectum.res)[1])
    #number of significantly enriched pathways
    paste0("The number of significantly enriched pathways is: ", num.pathways.rectum.sign <- dim(rectum.sign)[1])
    #exporting results to the file
    write.table(ewp.rectum.res, file=paste0("4-pathway_analysis/enrichResults_ORA_",disorder,"_rectum.tsv"),
                sep = "\t" ,quote = FALSE, row.names = FALSE)
    
  }
}

createHeatmap <- function(p_threshold_pathway, q_threshold_pathway){
  setwd(work_DIR)
  
  # Read files
  CD.ileum <- read.delim("4-pathway_analysis/enrichResults_ORA_CD_ileum.tsv",sep = "\t", header = TRUE)
  CD.rectum <- read.delim("4-pathway_analysis/enrichResults_ORA_CD_rectum.tsv", sep = "\t",header = TRUE)
  UC.ileum <- read.delim("4-pathway_analysis/enrichResults_ORA_UC_ileum.tsv",sep = "\t", header = TRUE)
  UC.rectum <- read.delim("4-pathway_analysis/enrichResults_ORA_UC_rectum.tsv", sep = "\t",header = TRUE)

  #we need to get pathways that has p.adjust value lower than 0.05 and qvalue<0.02
  #To prevent high false discovery rate (FDR) in multiple testing, q-values are also estimated for FDR control.
  CD.ileum.f <- CD.ileum[(CD.ileum$p.adjust<p_threshold_pathway)&(CD.ileum$qvalue<q_threshold_pathway),]
  CD.rectum.f <- CD.rectum[(CD.rectum$p.adjust<p_threshold_pathway)&(CD.rectum$qvalue<q_threshold_pathway),]
  UC.ileum.f <- UC.ileum[(UC.ileum$p.adjust<p_threshold_pathway)&(UC.ileum$qvalue<q_threshold_pathway),]
  UC.rectum.f <- UC.rectum[(UC.rectum$p.adjust<p_threshold_pathway)&(UC.rectum$qvalue<q_threshold_pathway),]
  
  #Filter out unused columns 
  CD.ileum.f  <- CD.ileum.f[,c(2,6)]
  CD.rectum.f <- CD.rectum.f[,c(2,6)]
  UC.ileum.f  <- UC.ileum.f[,c(2,6)]
  UC.rectum.f <- UC.rectum.f[,c(2,6)]
  #first 20 max value of p.adjust pathways for each comparison: Note that if a dataset has less then 20 sign. changed PWs, less rows need to be selected (e.g. adapt c(1:20))
  cd_ileum_temp <- CD.ileum.f
  cd_rectum_temp <- CD.rectum.f
  uc_ileum_temp <- UC.ileum.f
  uc_rectum_temp <- UC.rectum.f
  if (nrow(CD.ileum.f) > 20){
    cd_ileum_temp <- CD.ileum.f[c(1:20),]
  }
  if (nrow(CD.rectum.f) > 20){
    cd_rectum_temp <- CD.rectum.f[c(1:20),]
  }
  if (nrow(UC.ileum.f) > 20){
    uc_ileum_temp <- UC.ileum.f[c(1:20),]
  }
  if (nrow(UC.rectum.f) > 20){
    uc_rectum_temp <- UC.rectum.f[c(1:20),]
  }
  all.pathways.1 <- merge(cd_ileum_temp, cd_rectum_temp,by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
  #merge UC pathways
  all.pathways.2 <- merge(uc_ileum_temp, uc_rectum_temp, by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
  #merge all of them
  all.pathways <- merge(all.pathways.1 , all.pathways.2 ,by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
  colnames(all.pathways) <- c("Description","CD.ileum.p.adjust","CD.rectum.p.adjust","UC.ileum.p.adjust","UC.rectum.p.adjust")
  
  #replace NA values with the values from the whole list
  #### for CD ileum
  #find pathways which does not occur in the filtered enriched pathway list of cd.ileum (p.adjust<0.05 & qvalue<0.02 )
  #because we will not take into account not sig. enriched pathways, we will assign value of 1 for their p.adjust
  notExist.CDileum <- setdiff(all.pathways$Description,CD.ileum.f$Description)
  all.pathways[all.pathways$Description %in% notExist.CDileum,]$CD.ileum.p.adjust <- 1
  #the rest NA values correspond to the sig.enriched pathways but not in the first 20 list.
  #so we will replace NA values with the p.adjust values from the whole list
  NA.indices <- which(is.na(all.pathways$CD.ileum.p.adjust), arr.ind = TRUE)
  allIDs <- all.pathways[NA.indices,]$Description
  df <- CD.ileum.f[CD.ileum.f$Description %in% allIDs,]
  df <- df[order(df$Description),]
  all.pathways[all.pathways$Description %in% df$Description,]$CD.ileum.p.adjust <- df$p.adjust
  #### for CD rectum
  #find pathways which does not occur in the filtered enriched pathway list of cd rectum (p.adjust<0.05 & qvalue<0.02 )
  notExist.CDrectum<- setdiff(all.pathways$Description,CD.rectum.f$Description)
  all.pathways[all.pathways$Description %in% notExist.CDrectum,]$CD.rectum.p.adjust <- 1
  #replacing NA values with the values from the whole list
  NA.indices <- which(is.na(all.pathways$CD.rectum.p.adjust), arr.ind = TRUE)
  allIDs <- all.pathways[NA.indices,]$Description
  df <- CD.rectum.f[CD.rectum.f$Description %in% allIDs,]
  df <- df[order(df$Description),]
  all.pathways[all.pathways$Description %in% df$Description,]$CD.rectum.p.adjust <- df$p.adjust
  #### for UC ileum
  #find pathways which does not occur in the filtered enriched pathway list of UC ileum (p.adjust<0.05 & qvalue<0.02 )
  notExist.UCileum <- setdiff(all.pathways$Description,UC.ileum.f$Description)
  all.pathways[all.pathways$Description %in% notExist.UCileum,]$UC.ileum.p.adjust <- 1
  #### for UC rectum
  #find pathways which does not occur in the filtered enriched pathway list of UC rectum (p.adjust<0.05 & qvalue<0.02 )
  notExist.UCrectum<- setdiff(all.pathways$Description,UC.rectum.f$Description)
  all.pathways[all.pathways$Description %in% notExist.UCrectum,]$UC.rectum.p.adjust <- 1
  #replacing NA values with the values from the whole list
  NA.indices <- which(is.na(all.pathways$UC.rectum.p.adjust), arr.ind = TRUE)
  allIDs <- all.pathways[NA.indices,]$Description
  df <- UC.rectum.f[UC.rectum.f$Description %in% allIDs,]
  df <- df[order(df$Description),]
  all.pathways[all.pathways$Description %in% df$Description,]$UC.rectum.p.adjust <- df$p.adjust
  
  row.names(all.pathways) <- all.pathways$Description
  all.pathways  <- all.pathways[,2:5]
  colnames(all.pathways) <- c("CD Ileum","CD Rectum","UC Ileum","UC Rectum")
  
  #create output folder if not exist
  if(!dir.exists("5-create_heatmap")){dir.create("5-create_heatmap")}
  
  ## Select a size to visualize the heatmap with (options; large or small)
  size_heatmap <- "large"
  ##Print labels large for paper, small for notebook:
  fontsize_row_l = 30 
  if (size_heatmap == "large") {
    fontsize_col_l = 30 
    fontsize_l = 30
    width_l =2000 
    height_l =2000 
    name_heatmap_file <- "5-create_heatmap/heatmap_log10_large.png"
  }else if(size_heatmap == "small"){ 
    fontsize_row_l = 10 
    fontsize_col_l = 10 
    width_l =1500 
    height_l =1500 
    fontsize_l = 10
    name_heatmap_file <- "5-create_heatmap/heatmap_log10_small.png"
  }else{print("Size not Recognised")}
  #normally darker value represent higher values light color represent smaller values
  #when we use rev function higher ones are represented by light color
  colMain <- colorRampPalette(rev(brewer.pal(9, "Blues")))(30)
  my_heatmap <- pheatmap(as.matrix(log10(all.pathways)), scale = "none", color = colMain , 
                         legend = TRUE , legend_breaks = c(0, -5, -10, -15, min(log10(all.pathways))), 
                         main = "", 
                         legend_labels = c("-log10(adj. p-value) \n", " -5", " -10", " -15", ""),
                         cellwidth = 80, treeheight_row = 200, fontsize = fontsize_l, fontsize_row= fontsize_row_l, 
                         fontsize_col = fontsize_col_l, cluster_rows = TRUE, cluster_cols = FALSE)
  
  #save obtained heatmap
  save_pheatmap_png <- function(x, filename, width = width_l, height = height_l) {
    png(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  ##p-values are visalized on a log10 scale, to make them more discriminatory.
  save_pheatmap_png(my_heatmap, name_heatmap_file)
}

networkAnalysis <- function(){
  setwd(work_DIR)
  #read all DEG data
  CD.ileum <- read.delim("4-pathway_analysis/DEGs_CD_ileum.tsv",sep = "\t", header = TRUE)
  CD.rectum <- read.delim("4-pathway_analysis/DEGs_CD_rectum.tsv", sep = "\t",header = TRUE)
  UC.ileum <- read.delim("4-pathway_analysis/DEGs_UC_ileum.tsv",sep = "\t", header = TRUE)
  UC.rectum <- read.delim("4-pathway_analysis/DEGs_UC_rectum.tsv", sep = "\t",header = TRUE)

  #Listing all up and down regulated genes separately for CD:
  CD.up.ileum   <-unique(CD.ileum[CD.ileum$log2FC_ileum > 0,])
  colnames(CD.up.ileum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_CD", "pvalue_CD")
  CD.down.ileum <-unique(CD.ileum[CD.ileum$log2FC_ileum < 0,])
  colnames(CD.down.ileum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_CD", "pvalue_CD")
  CD.up.rectum   <-unique(CD.rectum[CD.rectum$log2FC_rectum > 0,])
  colnames(CD.up.rectum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_CD", "pvalue_CD")
  CD.down.rectum <-unique(CD.rectum[CD.rectum$log2FC_rectum < 0,])
  colnames(CD.down.rectum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_CD", "pvalue_CD")
  #Listing all up and down regulated genes separately for UC:
  UC.up.ileum   <-unique(UC.ileum[UC.ileum$log2FC_ileum > 0,])
  colnames(UC.up.ileum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_UC", "pvalue_UC")
  UC.down.ileum <-unique(UC.ileum[UC.ileum$log2FC_ileum < 0,])
  colnames(UC.down.ileum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_UC", "pvalue_UC")
  UC.up.rectum   <-unique(UC.rectum[UC.rectum$log2FC_rectum > 0,])
  colnames(UC.up.rectum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_UC", "pvalue_UC")
  UC.down.rectum <-unique(UC.rectum[UC.rectum$log2FC_rectum < 0,])
  colnames(UC.down.rectum) <- c ("HGNC_symbol", "ENTREZ", "log2FC_UC", "pvalue_UC")
  
  #OVERLAP ILEUM
  # overlap genes between CD down and UC down
  merged.ileum.downCDdownUC <- merge(x=CD.down.ileum, y=UC.down.ileum, by=c('ENTREZ', 'HGNC_symbol'), all.x=FALSE, all.y=FALSE)
  # overlap genes between CD up and UC down
  merged.ileum.upCDdownUC <- merge(x=CD.up.ileum,y=UC.down.ileum, by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  # overlap genes between CD up and UC up
  merged.ileum.upCDupUC <- merge(x=CD.up.ileum,y=UC.up.ileum,by=c('ENTREZ', 'HGNC_symbol'), all.x=FALSE, all.y=FALSE)
  # overlap genes between CD down and UC up
  merged.ileum.downCDupUC <- merge(x=CD.down.ileum,y=UC.up.ileum,by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  #merge all DEG with corresponding logFC for both diseases
  DEG.overlapped_ileum <- rbind(merged.ileum.downCDdownUC, merged.ileum.upCDdownUC, merged.ileum.upCDupUC, merged.ileum.downCDupUC)
  if(!dir.exists("6-network_analysis")) dir.create("6-network_analysis")
  write.table(DEG.overlapped_ileum ,"6-network_analysis/DEG.overlapped_ileum",row.names=FALSE,col.names = TRUE,quote= FALSE, sep = "\t")
  
  # OVERLAP RECTUM
  # overlap genes between CD down and UC down
  merged.rectum.downCDdownUC <- merge(x=CD.down.rectum,y=UC.down.rectum,by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  # overlap genes between CD up and UC down
  merged.rectum.upCDdownUC <- merge(x=CD.up.rectum,y=UC.down.rectum,by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  # overlap genes between CD up and UC up
  merged.rectum.upCDupUC <- merge(x=CD.up.rectum,y=UC.up.rectum,by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  # overlap genes between CD down and UC up
  merged.rectum.downCDupUC <- merge(x=CD.down.rectum,y=UC.up.rectum,by=c('ENTREZ', 'HGNC_symbol'),all.x=FALSE, all.y=FALSE)
  #merge all DEG with corresponding logFC for both diseases
  DEG.overlapped_rectum <- rbind(merged.rectum.downCDdownUC, merged.rectum.upCDdownUC, merged.rectum.upCDupUC, merged.rectum.downCDupUC)
  write.table(DEG.overlapped_rectum ,"6-network_analysis/DEG.overlapped_rectum",row.names=FALSE,col.names = TRUE,quote= FALSE, sep = "\t")
  
  # PPI NETWORK
  ## Select a location to analyse (options; ileum or rectum)
  for (location in c("ileum", "rectum")){
    if (location == "ileum"){deg <- DEG.overlapped_ileum}else if(location == "rectum"){deg <- DEG.overlapped_rectum}else{print("Location not recognized")}
    #check that cytoscape is connected
    cytoscapePing()
    #close session before starting --> this overwrites existing networks!
    closeSession(FALSE)
    networkName = paste0("PPI_network_",location)
    #create a PPI network using overlapped DEGs between CD and UC
    #first take the deg input list as a query
    x <- readr::format_csv(as.data.frame(deg$ENTREZ), col_names=F, escape = "double", eol =",")
    #then use the below function to convert a command string into a CyREST query URL, executes a GET request, 
    #and parses the result content into an R list object. Same as commandsGET
    commandsRun(paste0('string protein query cutoff=0.7', ' newNetName=',networkName, ' query=',x,' limit=0'))
    #get proteins (nodes) from the constructed network: #Query term: entrez IDs, display name: HGNC symbols.
    proteins <- RCy3::getTableColumns(columns=c("query term", "display name"))
    #get edges from the network
    ppis     <- RCy3::getTableColumns(table="edge", columns=c("name"))
    #split extracted edge information into source-target format
    ppis     <- data.frame(do.call('rbind', strsplit(as.character(ppis$name),' (pp) ',fixed=TRUE)))
    #merge obtained nodes and edges to get entrez IDs for each source genes 
    ppis.2   <- merge(ppis, proteins, by.x="X1", by.y="display name", all.x=T)
    #change column names
    colnames(ppis.2) <- c("s", "t", "source")
    #merge again to add entrez IDs of target genes 
    ppis.3   <- merge(ppis.2, proteins, by.x="t", by.y="display name", all.x=T)
    colnames(ppis.3)[4] <-"target"
    #ppi3 stores interaction between all proteins so add new column represeting type of interaction
    ppis.3$interaction <- "PPI"
    #add col names to protein
    colnames(proteins) <- c("id","label")
    proteins$type <- "protein"
    ###############get all pathways from WIKIPATHWAYS #################
    ##Work with local file (for publication), or new download:
    pathway_data <- "local" #Options: local, new
    if (pathway_data == "local") {
      wp.hs.gmt <-list.files(work_DIR, pattern="wikipathways", full.names=FALSE)
      paste0("Using local file, from: ", wp.hs.gmt )
    }else if(pathway_data == "new"){ 
      #below code should be performed first to handle the ssl certificate error while downloading pathways 
      options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
      #downloading latest pathway gmt files for human 
      wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
      paste0("Using new data, from: ", wp.hs.gmt)}else{print("Pathway data type not recognized")
      }
    #all wp and gene information stored in wp2gene object
    wp2gene   <- rWikiPathways::readPathwayGMT(wp.hs.gmt)
    #filter out  pathways that does not consist of any differentially expressed genes 
    wp2gene.filtered <- wp2gene [wp2gene$gene %in% deg$ENTREZ,]
    #change column names 
    colnames(wp2gene.filtered)[3] <- c("source")
    colnames(wp2gene.filtered)[5] <- c("target")
    #add new column for representing interaction type
    wp2gene.filtered$interaction <- "Pathway-Gene"
    #store only wp information 
    pwy.filtered <- unique( wp2gene [wp2gene$gene %in% deg$ENTREZ,c(1,3)])
    colnames(pwy.filtered) <- c("label", "id")
    pwy.filtered$type <- "pathway"
    colnames(pwy.filtered) <- c("label","id", "type")
    #get genes 
    genes <- unique(deg[,c(1,2), drop=FALSE])
    genes$type <- "gene"
    colnames(genes) <- c("id","label","type")
    genes$id <- as.character(genes$id)
    #genes and pathways are separate nodes and they need to be merged
    nodes.ppi <- dplyr::bind_rows(genes,pwy.filtered)
    rownames(nodes.ppi) <- NULL
    edges.ppi <- unique(dplyr::bind_rows(ppis.3[,c(3,4,5)], wp2gene.filtered[,c(3,5,6)]))
    rownames(edges.ppi) <- NULL
    
    networkName = paste0("PPI_Pathway_Network_",location)
    ###########Create PPI-pathway network###
    RCy3::createNetworkFromDataFrames(nodes= nodes.ppi, edges = edges.ppi, title=networkName, collection=location)
    RCy3::loadTableData(nodes.ppi, data.key.column = "label", table="node", table.key.column = "label")
    RCy3::loadTableData(deg, data.key.column = "ENTREZ", table.key.column = "id")
    ###########Visual style#################
    RCy3::copyVisualStyle("default","ppi")#Create a new visual style (ppi) by copying a specified style (default)
    RCy3::setNodeLabelMapping("label", style.name="ppi")
    RCy3::lockNodeDimensions(TRUE, style.name="ppi")#Set a boolean value to have node width and height fixed to a single size value.
    #threshold is set based of differential expressed gene criteria
    data.values<-c(-0.58,0,0.58) 
    #red-blue color schema chosen
    node.colors <- c(brewer.pal(length(data.values), "RdBu"))
    #nodes are split to show both log2fc values for both diseases 
    RCy3::setNodeCustomHeatMapChart(c("log2FC_CD","log2FC_UC"), slot = 2, style.name = "ppi", colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))
    RCy3::setVisualStyle("ppi")
    # Saving output
    outputName = paste0("6-network_analysis/PPI_Pathway_Network_", wp.hs.gmt, location,".png")
    png.file <- file.path(getwd(), outputName)
    exportImage(png.file,'PNG', zoom = 500)
    
    # CLUSTERING
    cytoscapePing()
    #Install the Clustermaker2 app (if not available already)
    if("Clustermaker2" %in% commandsHelp("")) print("Success: the Clustermaker2 app is installed") else print("Warning: Clustermaker2 app is not installed. Please install the Clustermaker2 app before proceeding.")
    if(!"Clustermaker2" %in% commandsHelp("")){
      installApp("Clustermaker2")
    }
    networkName = paste0 ("PPI_Pathway_Network_",location)
    #to get network name of the location 
    networkSuid = getNetworkSuid(networkName)
    setCurrentNetwork(network=networkSuid)
    #create cluster command
    clustermaker <- paste("cluster mcl createGroups=TRUE showUI=TRUE network=SUID:",networkSuid, sep="")
    #run the command in cytoscape
    res <- commandsGET(clustermaker)
    #total number of clusters 
    cat("Total number of clusters for",location, as.numeric(gsub("Clusters: ", "", res[1])))
    #change pathway node visualization
    pathways <- RCy3::selectNodes(nodes="pathway", by.col = "type")
    RCy3::setNodeColorBypass(node.names = pathways$nodes, "#D3D3D3")
    RCy3::setNodeBorderWidthBypass(node.names = pathways$nodes, 10)
    #export image
    outputName = paste0 ("6-network_analysis/PPI_Pathway_Network_",location, wp.hs.gmt,"_clustered",".png")
    png.file <- file.path(getwd(), outputName)
    exportImage(png.file,'PNG', zoom = 500)
    #save session
    cys.file <- file.path(getwd(), "6-network_analysis/PPI_Pathway_Network_",wp.hs.gmt, location,"_clustered",".cys")
  }
  
}


#############################################Metabolomics analysis functions ###################################################
#filtering metabolomics data 
filteringMets <- function(metaData,mbxData){
  
  #filter out by data type and week number
  metaDataMBX <- subset(metaData, metaData$data_type == "metabolomics" )
  #we need to have the samples which has same visit number
  metaDataMBX<- subset(metaDataMBX, metaDataMBX$visit_num == 4)
  print("Samples filtered data type and visit number")
  metaDataMBX <- metaDataMBX %>% dplyr::select(External.ID,Participant.ID,diagnosis)
  #rename columns of metaDataMBX
  colnames(metaDataMBX) <- c("ExternalID","ParticipantID","disease" )
  
  #delete not used columns
  mbxData = subset(mbxData, select = -c(1,2,3,4,7) )
  
  ### row (metabolite) filtering ###
  #delete metabolite or row if it has NA or empty value for hmdbID
  mbxData<- mbxData[!(is.na(mbxData$HMDB...Representative.ID.) | mbxData$HMDB...Representative.ID.=="") , ]
  #remove rows which has hmdb as "redundant ion"
  mbxData<- mbxData[!(mbxData$HMDB...Representative.ID.=="redundant ion") , ]
  #remove character (asterisk) in some hmdb column values
  mbxData$HMDB...Representative.ID.<- stringr::str_replace(mbxData$HMDB...Representative.ID., '\\*', '')
  #Update HMDB IDs to new data structure
  mbxData$HMDB...Representative.ID.<- stringr::str_replace(mbxData$HMDB...Representative.ID., 'HMDB', 'HMDB00')
  #back up original mbxdata
  mbxData.b <- mbxData
  
  ### modify mbxData based on sample names given in metaData file (created with the criteria visit_num=4 )###
  #filter out mbxData columns (samples) based metaDataMBX externalIDs
  names.use <- names(mbxData)[ names(mbxData) %in% metaDataMBX$ExternalID]
  #update mbx data with used names
  mbxData <- mbxData [ ,names.use]
  #order data based on col names
  mbxData <- mbxData[ , order(names(mbxData))]
  
  #order metadata based on externalID
  metaDataMBX <- metaDataMBX[order(metaDataMBX$ExternalID),]
  
  #add HMDBID and Compound Name column to the mbx data
  mbxData <- cbind(mbxData.b$HMDB...Representative.ID., mbxData.b$Metabolite,mbxData)
  colnames(mbxData)[1] <- "HMDB.ID"
  colnames(mbxData)[2] <- "Compound.Name"
  
  #add disease labels to the mbx data
  diseaseLabels <- metaDataMBX$disease
  ##Add two empty strings to macth with additional column data.
  diseaseLabels <- append(diseaseLabels, "NA",after = 0)
  diseaseLabels <- append(diseaseLabels, "NA",after = 0)
  mbxData <- rbind(diseaseLabels, mbxData)
  
  # split up the data for UC and CD, include the control data nonIBD
  #write only CD_healthy comparison
  mbxDataCD <- mbxData[ ,(mbxData[1, ] == "CD" | mbxData[1, ] == "nonIBD")]
  mbxDataCD <- cbind(mbxData[,1:2],mbxDataCD)
  colnames(mbxDataCD)[1]="HMBDB.ID"
  colnames(mbxDataCD)[2] <- "Compound.Name"
  
  #write only UC versus nonIBD comparison
  mbxDataUC <- mbxData[ ,(mbxData[1, ] == "UC" | mbxData[1, ] == "nonIBD")]
  #add hmdb id again
  mbxDataUC <- cbind(mbxData[,1:2],mbxDataUC)
  colnames(mbxDataUC)[1]="HMBDB.ID"
  colnames(mbxDataUC)[2] <- "Compound.Name"
  
  print("Samples and metabolites filtering process finished")  
 
  #Metabolites with >50% NA values will be filtered
  print("Filtering metabolites with NA values started")
  
  ################# for CD disease #########
  #Merge column headers: disorder_patientID
  names(mbxDataCD) <- paste(mbxDataCD [1, ], names(mbxDataCD), sep = "_")
  mbxDataCD <- mbxDataCD [-1,]

  columns <- ncol(mbxDataCD)
  rowsData <- nrow(mbxDataCD)
  removeLines <- rowSums(is.na(mbxDataCD[,3:columns]))
  fifty_percent <- floor((columns)/2)
  
  CD_MissingDataCounted <- cbind(mbxDataCD, removeLines)
  CD_NoMissingData <- subset(CD_MissingDataCounted, removeLines <= fifty_percent)
  #Remove last column for further processing.
  CD_NoMissingData <- subset(CD_NoMissingData, select=-c(removeLines))
  
  #Convert intensity data to numeric values                         
  CD_NoMissingData[, c(3:columns)] <- apply(CD_NoMissingData[, c(3:columns)],2, function(x) as.numeric(as.character(x)))
  write.table(CD_NoMissingData, "7-metabolite_data_preprocessing/filtered/mbxDataCD_nonIBD.csv", sep =",", row.names = FALSE)
  
  ############# for UC disease ################
  #Merge column headers: disorder_patientID
  names(mbxDataUC) <- paste(mbxDataUC [1, ], names(mbxDataUC), sep = "_")
  mbxDataUC <- mbxDataUC [-1,]
  
  columns <- ncol(mbxDataUC)
  rowsData <- nrow(mbxDataUC)
  removeLines <- rowSums(is.na(mbxDataUC[,3:columns])) #HERE WE DELETE A COLUMN WHICH SHOULD NOT BE LIKE THAT -it should be like-> mSet[,3:columns]
  fifty_percent <- floor((columns)/2)
  
  UC_MissingDataCounted <- cbind(mbxDataUC, removeLines)
  UC_NoMissingData <- subset(UC_MissingDataCounted, removeLines <= fifty_percent)
  #Remove last column for further processing.
  UC_NoMissingData <- subset(UC_NoMissingData, select=-c(removeLines))
  
  #Convert intensity data to numeric values                         
  UC_NoMissingData[, c(3:columns)] <- apply(UC_NoMissingData[, c(3:columns)],2, function(x) as.numeric(as.character(x)))
  write.table(UC_NoMissingData, "7-metabolite_data_preprocessing/filtered/mbxDataUC_nonIBD.csv", sep =",", row.names = FALSE)
  
  return(list((CD_NoMissingData),(UC_NoMissingData)))
}


#normalization metabolomics data 
normalizeMets <- function(selected = "log2 transformation"){

  #browser()
  splitted <- strsplit(selected, " ")
  transformation <- splitted [[1]][1]
  cat("Selected transformation ", transformation)
  
  #data will be read from filtered folder
  mSet_CD <- read.csv("7-metabolite_data_preprocessing/filtered/mbxDataCD_nonIBD.csv", na.strings=c("", "NA"))
  mSet_UC <- read.csv("7-metabolite_data_preprocessing/filtered/mbxDataUC_nonIBD.csv", na.strings=c("", "NA"))
  
  ############################# for CD disease #########################
  columns <- ncol(mSet_CD)
  if(transformation == "cube"){
    mSet_transformed <- cbind(mSet_CD[,c(1,2)], mSet_CD[,3:columns]^(1/3))
  }else if(transformation == "square"){
    mSet_transformed <- cbind(mSet_CD[,c(1,2)], mSet_CD[,3:columns]^(1/2))
  }else if(transformation == "log2"){
    mSet_transformed <- cbind(mSet_CD[,c(1,2)], log2(mSet_CD[,3:columns]))
  }else if(transformation == "log10"){
    mSet_transformed <- cbind(mSet_CD[,c(1,2)], log10(mSet_CD[,3:columns]))
  }else{print("Warning: name for transformation not recognized")}
  
  colnames(mSet_transformed)[1]="HMDB.ID"
  colnames(mSet_transformed)[2]="Compound.Name"
 
  ## Visualize the data after the transformation (for one sample to get an idea of suitability of transformation:
  #create histogram for original distribution for first column with data
  png(paste0("7-metabolite_data_preprocessing/normalized/CD_histogram_raw",".png"),width=1000,height=1000)
  hist(mSet_CD[,3], col='steelblue', main='Original')
  dev.off()
  cat("raw CD histogram created\n")
  
  #create histogram for log-transformed distribution 
  png(paste0("7-metabolite_data_preprocessing/normalized/CD_histogram_norm",".png"),width=1000,height=1000)
  hist(mSet_transformed[,3], col='steelblue', main='Original')
  dev.off()
  cat("norm CD histogram created\n")
  
  ########### Testing if the transformation creates a normally distributed dataset (alpha >= 0.05)
  ##Calculate all Shapiro values for raw and transformed data:
  mSet_NoMissingData_Shapiro <- lapply(mSet_CD[,3:columns], shapiro.test)
  mSet_transformed_Shapiro <- lapply(mSet_transformed[,3:columns], shapiro.test)
  
  #Obtain the p-values for raw and transformed data
  mSet_NoMissingData_Shapiro_pvalues <- do.call(rbind, mSet_NoMissingData_Shapiro)
  mSet_transformed_Shapiro_pvalues <- do.call(rbind, mSet_transformed_Shapiro)
  
  ## Count how often the p-value is above 0.05, to obtain an estimate of achieved normality due to transformation
  mSet_NoMissingData_Shapiro_pvalues_sum <- sum(mSet_NoMissingData_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)
  mSet_transformed_Shapiro_pvalues_sum <- sum(mSet_transformed_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)
  
  eighty_percent <- floor(((columns)/10)*8)
  
  CD_shapiro <- FALSE
  #Print relevant information:
  if(mSet_transformed_Shapiro_pvalues_sum[1] > eighty_percent )
  {
    paste0("Data after ", transformation ," transformation seems to follow a normal distribution for more then 80% of your data")
    CD_shapiro <- TRUE
  } 
  else{
    print("Advised to select a different data transformation procedure")
  }
  write.table(mSet_transformed, "7-metabolite_data_preprocessing/normalized/CD_norm_data.csv", sep =",", row.names = FALSE)
  
  
  ############################# for UC disease #########################
  columns <- ncol(mSet_UC)
  if(transformation == "cube"){
    mSet_transformed <- cbind(mSet_UC[,c(1,2)], mSet_UC[,3:columns]^(1/3))
  }else if(transformation == "square"){
    mSet_transformed <- cbind(mSet_UC[,c(1,2)], mSet_UC[,3:columns]^(1/2))
  }else if(transformation == "log2"){
    mSet_transformed <- cbind(mSet_UC[,c(1,2)], log2(mSet_UC[,3:columns]))
  }else if(transformation == "log10"){
    mSet_transformed <- cbind(mSet_UC[,c(1,2)], log10(mSet_UC[,3:columns]))
  }else{print("Warning: name for transformation not recognized")}
  
  colnames(mSet_transformed)[1]="HMDB.ID"
  colnames(mSet_transformed)[2]="Compound.Name"
  
  ## Visualize the data after the transformation (for one sample to get an idea of suitability of transformation:
  #create histogram for original distribution for first column with data
  png(paste0("7-metabolite_data_preprocessing/normalized/UC_histogram_raw",".png"),width=800,height=800)
  hist(mSet_UC[,3], col='steelblue', main='Original')
  dev.off()
  cat("raw UC histogram created\n")
  
  #create histogram for log-transformed distribution 
  png(paste0("7-metabolite_data_preprocessing/normalized/UC_histogram_norm",".png"),width=800,height=800)
  hist(mSet_transformed[,3], col='steelblue', main='Original')
  dev.off()
  cat("norm UC histogram created\n")
  
  ########### Testing if the transformation creates a normally distributed dataset (alpha >= 0.05)
  ##Calculate all Shapiro values for raw and transformed data:
  mSet_NoMissingData_Shapiro <- lapply(mSet_UC[,3:columns], shapiro.test)
  mSet_transformed_Shapiro <- lapply(mSet_transformed[,3:columns], shapiro.test)
  
  #Obtain the p-values for raw and transformed data
  mSet_NoMissingData_Shapiro_pvalues <- do.call(rbind, mSet_NoMissingData_Shapiro)
  mSet_transformed_Shapiro_pvalues <- do.call(rbind, mSet_transformed_Shapiro)
  
  ## Count how often the p-value is above 0.05, to obtain an estimate of achieved normality due to transformation
  mSet_NoMissingData_Shapiro_pvalues_sum <- sum(mSet_NoMissingData_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)
  mSet_transformed_Shapiro_pvalues_sum <- sum(mSet_transformed_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)
  
  eighty_percent <- floor(((columns)/10)*8)
  
  UC_shapiro <- FALSE
  #Print relevant information:
  if(mSet_transformed_Shapiro_pvalues_sum[1] > eighty_percent )
  {
    paste0("Data after ", transformation ," transformation seems to follow a normal distribution for more then 80% of your data")
    UC_shapiro <- TRUE
  } 
  else{
    print("Advised to select a different data transformation procedure")
  }
  write.table(mSet_transformed, "7-metabolite_data_preprocessing/normalized/UC_norm_data.csv", sep =",", row.names = FALSE)
  
  return (list(CD_shapiro, UC_shapiro))
}

####################################################STATISTICAL ANALYSIS FUNCTIONS  #########################################

statAnalysisMets <- function (mSet_transformed,disorder,transformation, FC, pvalue){
  
 # browser()
  
  #take only first token of the transformation 
  transformation <- strsplit(transformation, " ")[[1]][1]
  
  if (disorder == "CD") {
    mSet_FINAL <- mSet_transformed[ , order(names(mSet_transformed))]}
  else{
    mSet_FINAL <- mSet_transformed[ , order(names(mSet_transformed), decreasing=TRUE)]}
 
  ##Move Name column back to column 2
  columnNumber <- which(colnames(mSet_FINAL)=="Compound.Name")
 # mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL)-1)]
  mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL))]
  ## WHILE APPLYING THIS STEP ONE COLOUMN WAS DELETED ->  FOR CD "nonIBD_MSM9VZMM" is deleted
  #we need to change it as mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL))
  mSet_FINAL <- mSet_FINAL [,-(columnNumber+1)]
  
  ##Move HMDB column back to the start
  columnNumber <- which(colnames(mSet_FINAL)=="HMDB.ID")
  #mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL)-1)]
  mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL))]
  mSet_FINAL <- mSet_FINAL [,-(columnNumber+1)]
  
  #Find relevant columns per group.
  library(stringr)
  columns_disorders <- sum(str_count(colnames(mSet_FINAL), disorder))
  end_Disorders <- (columns_disorders+2) ##add two, to compensate for HMDB.ID and Compound.Names column.
  
  ##calculate logFC for 2 groups (CD vs control, UC vs control), ignoring missing values (NAs) when calculating the mean.  
  disease = apply(mSet_FINAL[,3:end_Disorders], 1, mean, na.rm=TRUE)
  control_IBD = apply(mSet_FINAL[,(end_Disorders+1):ncol(mSet_FINAL)], 1, mean, na.rm=TRUE)
  
  #because the metabolomics data is already log2 transformed, we need to take the difference between the means (iso dividing the means over one another), since log2 Fold Change or log2 Ratio == log2(condition / control). Note: if the transformation step applied is cube_root or square_root, one needs to divide control over disease for this step!
  if(transformation == "log2" | transformation == "log10"){
    foldchange_disorder <-  disease - control_IBD
  }else{
    foldchange_disorder <- (disease /control_IBD )}
  
  ##ADD HMDB column at start, add fold change columns.
  mSet_AnalysisReady <- cbind(mSet_FINAL$HMDB.ID, mSet_FINAL$Compound.Name, foldchange_disorder)
  ##Rename first column and second column
  colnames(mSet_AnalysisReady)[1] <- "HMDB_ID"
  colnames(mSet_AnalysisReady)[2] <- "Compound_Name"
  
  ##Calculate p-value for two groups based on t-test (comparing control to disease).
  ##general function to store p-values for multiple rows:
  ttest_mSet <- function(df, grp1, grp2) {
    x = df[grp1]
    y = df[grp2]
    x = as.numeric(x)
    y = as.numeric(y)  
    results = t.test(x, y)
    results$p.value
  }
  p_values_disorder <- apply(mSet_FINAL, 1, ttest_mSet, grp1 = c(3:end_Disorders), grp2 = c((end_Disorders+1):ncol(mSet_FINAL)))
  
  ##Add p_values column to analysis dataset:
  mSet_AnalysisReady <- cbind(mSet_AnalysisReady, p_values_disorder)
  
  #Convert logFC and p-values columns to numeric values            
  mSet_AnalysisReady <- as.data.frame(mSet_AnalysisReady)
  mSet_AnalysisReady[ , c(3,4)] <- apply(mSet_AnalysisReady[ , c(3,4)], 2, function(x) as.numeric(as.character(x)))

  ## Volcano plot
  #The next step will visualize relevant metabolites in a Volcano plot, with thresholds defined for the log2FCs and p-values.
 # Since we are dealing with non-targeted LC_MS data, we are applying stringent criteria here (p-value below 0.05, and |log2FC| >= 1.5)
 
  ##Inspired by: https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
  if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2")}
  library('ggplot2')
  
  ##Define the thresholds for log2 (Fold Change) and p-values
  #For cut-off value uncertainties, see https://doi.org/10.1039/C6AN01342B .
  #log2FC_min <- -log2(FC)
  log2FC_min <- as.numeric(format(round(-log2(FC), 2), nsmall = 2))
  #log2FC_max <- log2(FC)
  log2FC_max <- as.numeric(format(round(log2(FC), 2), nsmall = 2))
  p_value_threshold <- pvalue
  
  ##Create column with HMDB_IDs, only if the data is relevant
  mSet_AnalysisReady$relevant_labels <- mSet_AnalysisReady$HMDB_ID
  mSet_AnalysisReady$relevant_labels[!((mSet_AnalysisReady$foldchange_disorder <= log2FC_min | mSet_AnalysisReady$foldchange_disorder >= log2FC_max) &  mSet_AnalysisReady$p_values_disorder <= p_value_threshold)] <- NA
  
  ##Duplication issues:
  if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}
  library(dplyr)
  
  ##Check for duplicate HMDB.IDs, however with different compound name (example: HMDB0011130, called LPE, LPE-A (both relevant), and LPE-B (not relevant) & HMDB0000252, called sphingosine isomer 1, 2 and 3.)
  mSet_AnalysisReady_Duplicates <- mSet_AnalysisReady %>% group_by(HMDB_ID) %>% distinct(Compound_Name, .keep_all = TRUE)
  #Remove all rows with relevant data, which will be added later by their average in case of duplicate HMDB.IDs:
  mSet_AnalysisReady_Duplicates <- subset(mSet_AnalysisReady_Duplicates, is.na(mSet_AnalysisReady_Duplicates$relevant_labels))
  mSet_AnalysisReady_Duplicates <- subset(mSet_AnalysisReady_Duplicates[,1:4])
  
  ##Check for duplicates HMDB.IDs, calculate the average for these if both significant (example: HMDB0010384).
  mSet_AnalysisReady_FC <- mSet_AnalysisReady %>% group_by(HMDB_ID, Compound_Name) %>% filter(!is.na(relevant_labels)) %>% summarize(foldchange_disorder=mean(foldchange_disorder)) 
  mSet_AnalysisReady_p <- mSet_AnalysisReady %>% group_by(HMDB_ID, Compound_Name)  %>% filter(!is.na(relevant_labels)) %>% summarize(p_values_disorder=mean(p_values_disorder))
  
  #Merge FC and p-value data together based on HMDB.ID
  mSet_AnalysisReady_FCandp <-
    left_join(mSet_AnalysisReady_FC, mSet_AnalysisReady_p, by=c("HMDB_ID", "Compound_Name")) %>%
    rowwise()  
  colnames(mSet_AnalysisReady_FCandp)[3] <- "foldchange_disorder"
  colnames(mSet_AnalysisReady_FCandp)[4] <- "p_values_disorder"
  
  ##Combine datasets again, before continuing to next steps:
  mSet_AnalysisFinal <- rbind(mSet_AnalysisReady_Duplicates, mSet_AnalysisReady_FCandp)
  
  ##Create column with HMDB_IDs again, only if the data is relevant, for visualization in Volcano Plot
  mSet_AnalysisFinal$relevant_ids <- mSet_AnalysisFinal$HMDB_ID
  mSet_AnalysisFinal$relevant_ids[!((mSet_AnalysisFinal$foldchange_disorder <= log2FC_min | mSet_AnalysisFinal$foldchange_disorder >= log2FC_max) &  mSet_AnalysisFinal$p_values_disorder <= p_value_threshold)] <- NA
  
  ##Create another column with Compound names, only if the data is relevant, for visualization in Volcano Plot
  mSet_AnalysisFinal$relevant_labels <- mSet_AnalysisFinal$Compound_Name
  mSet_AnalysisFinal$relevant_labels[!((mSet_AnalysisFinal$foldchange_disorder <= log2FC_min | mSet_AnalysisFinal$foldchange_disorder >= log2FC_max) &  mSet_AnalysisFinal$p_values_disorder <= p_value_threshold)] <- NA
  
  ##Select visualization of HMDB-IDs(relevant_ids), or chemical compound names(relevant_labels)
  selectViz <- "relevant_labels"
  
  # Library needed to space out the labels
  if(!"ggrepel" %in% installed.packages()){install.packages("ggrepel")}
  library(ggrepel)
  
  if(selectViz == "relevant_labels"){
    ##volcanoPlot_Disorder 
    volcanoPlot_disorder <- ggplot(data=mSet_AnalysisFinal, aes(x=foldchange_disorder, y=-log10(p_values_disorder), label=relevant_labels)) + geom_point() + theme_minimal() + geom_text_repel()
  }else if(selectViz == "relevant_ids"){
    ##volcanoPlot_Disorder 
    volcanoPlot_disorder <- ggplot(data=mSet_AnalysisFinal, aes(x=foldchange_disorder, y=-log10(p_values_disorder), label=relevant_ids)) + geom_point() + theme_minimal() + geom_text_repel()
  }else{print("Column name not recognized for label visualization.")}
  
  ## Add vertical lines for FoldChange and P-value thresholds:
  volcanoPlot_disorder <- volcanoPlot_disorder + geom_vline(xintercept=c(log2FC_min, log2FC_max), col="blue") +
    geom_hline(yintercept=-log10(p_value_threshold), col="red") + theme(plot.background = element_rect(fill = "white"))
  
  if (disorder == "UC") {disorderName <- "Ulcerative Colitis (UC)"}else{disorderName <- "Crohn's disease (CD)"}
  
  titleVolcano <- paste0("Volcano plot of ", transformation, " transformed data for ", disorderName )
  verticalAxisTitle <- paste0(transformation, " Fold Change, ", disorderName, " versus control ")
  
  ## Add title and update axis labels:
  volcanoPlot_disorder <- volcanoPlot_disorder + ggtitle(titleVolcano) + labs(y = "-log10(p-value)", x = verticalAxisTitle)
  
  ## Export the data and Volcano Plot:
  
  ##Save the data file
  nameDataFile <- paste0("8-significantly_changed_metabolites_analysis/mbxData_", disorder ,".csv")
  write.table(mSet_AnalysisFinal, nameDataFile, sep =",", row.names = FALSE)
  
  if(!"svglite" %in% installed.packages()){install.packages("svglite")}
  library('svglite')
  
  ##Save the Volcano plot:
  imageType <- "png" ##Options are: svg, png, eps, ps, tex, pdf, jpeg, tiff, png, bmp, svg or wmf
  #nameVolcano <- paste0("8-significantly_changed_metabolites_analysis/", disorder, "_", selectViz, "_VolcanoPlot_absLogFC_", log2FC_max, "_pValue_", p_value_threshold, ".", imageType)
  nameVolcano <- paste0("8-significantly_changed_metabolites_analysis/", disorder, "_", selectViz, "_VolcanoPlot.", imageType)
  ggsave(nameVolcano)
  
  
  
}

#identifier mapping for metabolomics data
mappingMets <- function (){
 
  #to set timeout option for downloading data 
  options(timeout=10000)

  #Obtain data from step 8
  mSet_CD <- read.csv("8-significantly_changed_metabolites_analysis/mbxData_CD.csv", na.strings=c("", "NA"))
  mSet_UC <- read.csv("8-significantly_changed_metabolites_analysis/mbxData_UC.csv", na.strings=c("", "NA"))
  #filter out unused columns
  mSet_CD <- mSet_CD [,c(1:4)]
  mSet_UC <- mSet_UC [,c(1:4)]
  
  ##Retain all metabolite IDs from CD and UC for identifier mapping:
  mSet_total <- unique(rbind(mSet_CD[,c(1,2)], mSet_UC[,c(1,2)]))
  
  ##Set the working directory to download the Metabolite mapping file
  location <- "data/metabolites.bridge"
  checkfile <- paste0(getwd(), '/' ,location)
  
  ##Download the Metabolite mapping file (if it doesn't exist locally yet):
  if (!file.exists(checkfile)) {
    download.file("https://figshare.com/ndownloader/files/26001794",location) #this code should be changed since the full size of the file can not 
    #be downloaded
  }
  
  mapper <- BridgeDbR ::loadDatabase(checkfile)
  
  ## Obtain the System codes for the databases HMDB (source database of dataset) and ChEBI (intended output database)
  code <- getOrganismCode("Homo sapiens")
  code_mappingFrom <- getSystemCode("HMDB")
  code_mappingTo   <- getSystemCode("ChEBI")
  
  ## Create a data frame with the mappings and the correct SystemCode
  input = data.frame(
    source = rep(code_mappingFrom, length(mSet_total$HMDB_ID)),
    identifier = mSet_total$HMDB_ID)
  #Obtain all mappings from HMDB to ChEBI
  MultiMappings = BridgeDbR::maps(mapper, input, code_mappingTo)
  #remove all rows in the mapped data which do not include the prefix "CHEBI"
  MultiMappings <- MultiMappings %>% filter(grepl("CHEBI",mapping, fixed = TRUE))
  #filter out double identifiers because there are one-to-many relationship between hmdb and chebi IDs
  MultiMappings <- MultiMappings %>% distinct (MultiMappings$identifier, .keep_all = TRUE)
  MultiMappings <- MultiMappings [,c(2,4)]
  
  merged.data_CD<- merge(MultiMappings, mSet_CD,by.x="identifier", by.y="HMDB_ID",sort = TRUE, all.x = TRUE, all.y = TRUE)
  merged.data_UC<- merge(MultiMappings, mSet_UC,by.x="identifier", by.y="HMDB_ID",sort = TRUE, all.x = TRUE, all.y = TRUE)
  #filter out metabolites that has NA value for CHEBI
  merged.data_CD<- merged.data_CD %>% tidyr::drop_na(mapping)
  merged.data_UC<- merged.data_UC %>% tidyr::drop_na(mapping)
  #filter out metabolites that has NA value for foldchange_disorder
  merged.data_CD<- merged.data_CD %>% tidyr::drop_na(foldchange_disorder)
  merged.data_UC<- merged.data_UC %>% tidyr::drop_na(foldchange_disorder)
  #change column names
  colnames(merged.data_CD) <- c("HMDBID","CHEBI", "label", "log2FC_met", "pvalue_met")
  colnames(merged.data_UC) <- c("HMDBID","CHEBI", "label", "log2FC_met", "pvalue_met")

  
  ## Export the mapped metabolomics data:
  if(!dir.exists("9-metabolite_identifier_mapping"))
    dir.create("9-metabolite_identifier_mapping")
  ##Save the data file
  write.csv(merged.data_CD, '9-metabolite_identifier_mapping/mbx_mapped_data_CD.csv', row.names = FALSE)
  write.csv(merged.data_UC, '9-metabolite_identifier_mapping/mbx_mapped_data_UC.csv', row.names = FALSE)

}

#pathway analysis for metabolomics data 
pathwayAnalysisMets <- function (mSet,disorder){
 
  if (disorder == "CD") {
    print("Selected disorder is Crohn's disease")}
  else if(disorder == "UC"){ 
      print("Selected disorder is Ulcerative Colitis")}
  else{print("Disorder not Recognised")}
  
  
  if(!"SPARQL" %in% installed.packages()){
    install.packages("SPARQL")
  }
  library(SPARQL)
  ##Connect to Endpoint WikiPathways
  endpointwp <- "https://sparql.wikipathways.org/sparql"
  ## 1. Query metadata:
  queryMetadata <-
    "SELECT DISTINCT ?dataset (str(?titleLit) as ?title) ?date ?license 
      WHERE {
         ?dataset a void:Dataset ;
         dcterms:title ?titleLit ;
         dcterms:license ?license ;
         pav:createdOn ?date .
       }"
  #below code should be performed first to handle the ssl certificate error
  options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
  resultsMetadata <- SPARQL(endpointwp,queryMetadata,curl_args=list(useragent=R.version.string))
  showresultsMetadata <- resultsMetadata$results
  remove(queryMetadata, resultsMetadata)
   
  
  ## Create a list of HMDB IDs according to filtering criteria from step 8.
  list_Relevant_HMDB_IDs <- list(mSet$relevant_ids)
  vector_HMDB <- unlist(list_Relevant_HMDB_IDs) #convert list to array, for traversing the data to a SPARQL query later on
  vector_HMDB <- vector_HMDB[!is.na(vector_HMDB)]
  ##Add the HMDb prefix IRI in front of all IDs.
  query_HMDBs <- paste("ch:", vector_HMDB, sep="")
  ##Merge the individual entries in the vector into one string, separated by a space
  string_HMDB <- paste(c(query_HMDBs), collapse=' ' )
  
  ##TODO: add column with nr. of Metabolites in PW (to calculate PW impact)
  
  #For now, filter out Reactome PWs due to visualization issues in Cytoscape.
  item1 = "PREFIX ch: <https://identifiers.org/hmdb/>
PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (count(distinct ?hmdbMetabolite) AS ?HMDBsInPWs) (count(distinct ?metaboliteDatanode) AS ?TotalMetabolitesinPW) (GROUP_CONCAT(DISTINCT fn:substring(?hgnc,37);separator=' ') AS ?Proteins) (count(distinct ?hgnc) AS ?ProteinsInPWs) (GROUP_CONCAT(DISTINCT fn:substring(?hmdbMetabolite,30);separator=' ') AS ?includedHMDBs)
where {
VALUES ?hmdbMetabolite {"
  item2 = "}
 
 ?metaboliteDatanode	a wp:Metabolite ;
                        dcterms:isPartOf ?pathwayRes .
 
 ?datanode	a wp:Metabolite ;          
           	wp:bdbHmdb  ?hmdbMetabolite ;
    		dcterms:isPartOf ?pathwayRes .
 ?pathwayRes a wp:Pathway ;
             wp:organismName 'Homo sapiens' ; 
    		dcterms:identifier ?wpid ;
    		dc:title ?title .
    		
 ?datanode2 wp:bdbHgncSymbol ?hgnc ;
    		dcterms:isPartOf ?pathwayRes .
    		
  #?pathwayRes wp:ontologyTag cur:Reactome_Approved . 
  ?pathwayRes wp:ontologyTag cur:AnalysisCollection .   		
}
ORDER BY DESC(?HMDBsInPWs)"
query_CombinePWs <- paste(item1,string_HMDB,item2)
remove(item1, item2)

results_CombinePWs <- SPARQL(endpointwp,query_CombinePWs,curl_args=list(useragent=R.version.string))
showresults_CombinePWs <- results_CombinePWs$results
remove(query_CombinePWs,results_CombinePWs)
#Print table with first 5 relevant pathways (if less than 5 are found, print only those)
if(nrow(showresults_CombinePWs) < 5){
  print(showresults_CombinePWs[1:nrow(showresults_CombinePWs),c(2:5)])
}else{print(showresults_CombinePWs[1:5,c(2:5)])}

remove(cleaned_string_HMDB, list_Relevant_HMDB_IDs, query_HMDBs)
  
## Calculate the ORA score for each pathway, using the Fishers exact test.
#Create a dataframe to store the required numbers in.
Contingency_table <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("WP.ID", "x", "m", "n", "k"))))
counter = 1
for (i in 1:nrow(showresults_CombinePWs)) {
  Contingency_table[counter,1] <- (showresults_CombinePWs[i,2]) #WP.ID
  Contingency_table[counter,2] <- (showresults_CombinePWs[i,4]) ##x <- (number4) #Total differentially changed metabolites, also in a PW. (HMDBsInPWs)
  Contingency_table[counter,3] <- (showresults_CombinePWs[i,5]) ##m <- (number) #Total Metabolites in PW (TotalMetabolitesinPW)
  Contingency_table[counter,4] <- (length(unique(mSet[,1])) - showresults_CombinePWs[i,4]) ##n <- (number2) #Total Metabolites measured not in PW (DISTINCT all_HMDB - HMDBsInPWs)
  Contingency_table[counter,5] <- length(unique(vector_HMDB)) ##k <- (number3) #Total differentially changed metabolites. (DISTINCT vector_HMDB)
  
  counter <- counter + 1
}

# Calculate hypergeometric density p-value for all pathways.
i <- 1:nrow(Contingency_table)
probabilities <- dhyper(Contingency_table[i,2], Contingency_table[i,3], Contingency_table[i,4], Contingency_table[i,5], log = FALSE)

pathwayAnalysis_results <- cbind(showresults_CombinePWs[, c(2:4)], probabilities, showresults_CombinePWs[, c(6,7)])
colnames(pathwayAnalysis_results)[5] <- "HGNCs"
colnames(pathwayAnalysis_results)[6] <- "ProteinsInPWs"

##Sort PW results based on 1. highest amount of #HMDBs in PW, 2. lowest p-values,  3. highest amouny of proteins in PW (which might be relevant for transcriptomics analysis later)
pathwayAnalysis_results_sorted <- pathwayAnalysis_results[  with(pathwayAnalysis_results, order(-HMDBsInPWs, probabilities, -ProteinsInPWs)),]

print(pathwayAnalysis_results_sorted[1:5,])

if(!dir.exists("10-metabolite_pathway_analysis"))
  dir.create("10-metabolite_pathway_analysis")

## Export the pathway data:
nameDataFile <- paste0("10-metabolite_pathway_analysis/mbxPWdata_", disorder ,".csv")
write.table(pathwayAnalysis_results_sorted, nameDataFile, sep =",", row.names = FALSE)

}


##########################
## saveStatOutputDESeq2 ##
##########################

saveStatOutputDESeq2 <- function(design,dds,annotation=NULL,na.rm=TRUE,includeIntercept=FALSE,createHists=TRUE,postfix="") {
  files <- NULL
  for (i in 1:length(colnames(design))) {
    if(includeIntercept | length(grep("intercept",colnames(design)[i],ignore.case=TRUE))==0) {
      cat(paste("--[[ Saving table for coefficient ", colnames(design)[i], " ]]--\n", sep="\t"))
      res <- results(dds, contrast=design[,i], format="DataFrame")
      toptab <- as.data.frame(res)
      if(!is.null(annotation)) {
        toptab <- cbind(toptab,annotation[match(rownames(toptab),rownames(annotation)),])
      }
      fc <- as.matrix(2^toptab[,"log2FoldChange"])
      fc[(toptab[, "log2FoldChange"] < 0)&(!is.na(toptab[, "log2FoldChange"]))] <- -1/fc[(toptab[, "log2FoldChange"] < 0)&(!is.na(toptab[, "log2FoldChange"]))]
      colnames(fc) <- "FoldChange"
      m.col <- grep("log2FoldChange",colnames(toptab))
      toptab <- cbind(toptab[,1:m.col,drop=FALSE],fc,toptab[,(m.col+1):dim(toptab)[2],drop=FALSE])
      if(na.rm) {
        cat("----[[ ",sum(is.na(toptab[,"pvalue"]))," probes with NA estimates removed ]]\n")
        toptab <- toptab[!is.na(toptab[,"pvalue"]),]
      }
      filename <- paste("table_",colnames(design)[i],ifelse(postfix!="","_",""),postfix,".tab",sep="")
      write.table(toptab,file=filename,sep="\t",col.names=NA)
      files <- c(files, filename)
      if (createHists) {
        #also save p value histogram
        png(paste("pvalue_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
        hist(toptab[,"pvalue"],main=paste("p-value histogram for",colnames(design)[i]),xlab="p-values",col="blue",breaks=120,cex.axis=1.2,cex.lab=1.2)
        dev.off()
        #and of adapted FC
        png(paste("FC_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
        hist(toptab[,"FoldChange"],main=paste("adapted fold change histogram for",colnames(design)[i]),xlab="adapted fold changes",col="green3",breaks=120,cex.axis=1.2,cex.lab=1.2)
        dev.off()
      }
    }
  }
  return(files)
}

#####################################
###Create significant pvalue table###
#####################################

createPvalTab <- function(files,postfix="",pvaluelist=c(0.001,0.01,0.05,0.1),
                          adjpvaluelist=0.05,foldchangelist=c(1.1,1.2,1.5,2),
                          namePVal="P.Value",nameAdjPVal="adj.P.Val",
                          nameFC="Fold.Change",nameLogFC="logFC",html=FALSE, FC_threshold = 1) { ###### log2FC added as a user-variable

  if(is.null(files)) stop("vector of one or more file names has to be provided")
  
  if(html) library(R2HTML)
  
  #extract data from the comparison tables
  for(i in 1:length(files)){
    #read file
    tab <- read.delim(files[i],header=TRUE,as.is=TRUE)
    
    #create headers
    if(i==1) {
      numberOfGenes <- dim(tab)[1]
      ### P-values
      #expected
      expPval <- as.vector(apply(as.data.frame(pvaluelist),1,function(x) { 
        c(y <- ceiling(numberOfGenes*x),ceiling(0.5 * y),ceiling(0.5 * y))}))
      
      P.Values <- c("Comparisons",paste("pVal <",rep(pvaluelist,each=3),c("total","up","down")),"NAs")
      Expected <- c("Number expected",expPval,0)
      P.Values <- cbind(P.Values,Expected)
      
      ### Adjusted p-Values
      Adj.pValues <- c("Comparisons",paste("adj.pVal <",rep(adjpvaluelist,each=3),c("total","up","down")),"NAs")
      
      ### Fold changes
      Fold.Changes <- c("Comparisons",paste(rep(c("|FC| >=","FC >=","FC <="),3),rep(foldchangelist,each=3),c("total","up","down")),"NAs")
    }
    
    #remove extension from file name for display in table
    nameNoExt <- paste(strsplit(files[i],".",fixed=TRUE)[[1]][-length(strsplit(files[1],".",fixed=TRUE)[[1]])],collapse=".")
    
    #log2 of selected FC value   
    log2FC <- log2(FC_threshold)
    
    ### P-values
    rowsP <- NULL
    for (count in 1:length(pvaluelist)){
      noOfPval <- up <- down <- 0
      #compare p-value with entered p-values 
      noOfPval <- sum(tab[,namePVal] < pvaluelist[count],na.rm=TRUE)
      up <- sum((tab[,namePVal] < pvaluelist[count]) & (tab[,nameLogFC] >= log2FC),na.rm=TRUE) ### this should be user variable
      down <- sum((tab[,namePVal] < pvaluelist[count]) & (tab[,nameLogFC] <= (-log2FC)),na.rm=TRUE)
      rowsP <- c(rowsP,noOfPval,up,down)
    }
    P.Values<-cbind(P.Values,c(nameNoExt,rowsP,sum(is.na(tab[,namePVal]))))
    
    ### Adjusted p-Values
    rowsAP <- NULL
    for (count in 1:length(adjpvaluelist)){
      adjNoOfPval <- adjUp <- adjDown <- 0
      #compare adjusted p-values with entered adjusted p-values 
      adjNoOfPval <- sum(tab[,nameAdjPVal] < adjpvaluelist[count],na.rm=TRUE)	
      adjUp <- sum((tab[,nameAdjPVal] < adjpvaluelist[count]) & (tab[,nameLogFC] >= log2FC),na.rm=TRUE)
      adjDown <- sum((tab[,nameAdjPVal] < adjpvaluelist[count]) & (tab[,nameLogFC] <= (-log2FC) ),na.rm=TRUE)
      rowsAP <- c(rowsAP,adjNoOfPval,adjUp,adjDown)
    }
    Adj.pValues<-cbind(Adj.pValues,c(nameNoExt,rowsAP,sum(is.na(tab[,nameAdjPVal]))))
    
    ### Fold changes
    rowsFC <- NULL
    for (count in 1:length(foldchangelist)){
      FCTot <- FCUp <- FCDown <- 0
      #compare FC with entered FC 
      FCTot <- sum(abs(tab[,nameFC]) >= foldchangelist[count] & (!is.na(tab[,nameFC])))
      FCUp <- length(tab[(tab[,nameFC] >= foldchangelist[count]) & (!is.na(tab[,nameFC])),FC_threshold])
      FCDown <- length(tab[(tab[,nameFC] <= (-foldchangelist[count])) & (!is.na(tab[,nameFC])),FC_threshold])
      rowsFC <- c(rowsFC,FCTot,FCUp,FCDown)
    }
    Fold.Changes<-cbind(Fold.Changes,c(nameNoExt,rowsFC,sum(is.na(tab[,nameFC]))))
    
  }
  
  colnames(P.Values) <- NULL
  colnames(Adj.pValues) <- NULL
  colnames(Fold.Changes) <- NULL
  
  #write table to tab delimited text file 
  filename <- paste("Summary_tables",ifelse(postfix!="","_",""),postfix,".tab",sep="")
  cat (paste("--[[ Saving",filename,"]]--\n"))
  write.table(t(cbind(c("Total number of measurements",numberOfGenes),rep("",2))),file=filename,sep="\t",
              row.names=FALSE,col.names=FALSE, quote=FALSE)
  write.table(t(P.Values),file=filename,sep="\t",
              row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  write.table(t(cbind(rep("",dim(Adj.pValues)[1]),Adj.pValues)),file=filename,sep="\t",
              row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  write.table(t(cbind(rep("",dim(Fold.Changes)[1]),Fold.Changes)),file = filename,sep="\t",
              row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  
  #if requested, write table to HTML file
  if(html) {
    filenameHTML <- paste("Summary_tables",ifelse(postfix!="","_",""),postfix,".html",sep="")  
    cat (paste("--[[ Saving",filenameHTML,"]]--\n"))
    HTML(t(c("Total number of measurements",numberOfGenes)),file=filenameHTML,
         row.names=TRUE,append=FALSE)
    HTML(t(P.Values),file=filenameHTML,row.names=TRUE)
    HTML(t(Adj.pValues),file=filenameHTML,row.names=TRUE)
    HTML(t(Fold.Changes),file=filenameHTML,row.names=TRUE)
  }
  
}

#############################################Multi-omics analysis functions ###################################################

pathwaySelection <- function(p_threshold_multi_trans,
                             q_threshold_multi_trans,
                             p_threshold_multi_mets,
                             nProteinsPathway,
                             nMetsPathway){
  filelocation_t <- paste0(work_DIR, "/4-pathway_analysis/")
  #Obtain data from step 4 (transcript PWs)
  tPWs_CD_ileum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_CD_ileum.tsv'), sep = "\t", header = TRUE)
  tPWs_CD_rectum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_CD_rectum.tsv'), sep = "\t",header = TRUE)
  tPWs_UC_ileum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_UC_ileum.tsv'), sep = "\t", header = TRUE)
  tPWs_UC_rectum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_UC_rectum.tsv'), sep = "\t",header = TRUE)
  #Set location to download data for metabolomics pathway analysis:
  filelocation_m <- paste0(work_DIR, "/10-metabolite_pathway_analysis/")
  #Obtain data from step 9 (metabolite PWs)
  mPWs_CD <- read.delim(paste0(filelocation_m, 'mbxPWdata_CD.csv'), sep = ",", na.strings=c("", "NA"))
  mPWs_UC <- read.delim(paste0(filelocation_m, 'mbxPWdata_UC.csv'), sep = ",", na.strings=c("", "NA"))
  
  # Select pathways
  tPWs_CD_ileum_sign <- tPWs_CD_ileum[(tPWs_CD_ileum$p.adjust<p_threshold_multi_trans)&(tPWs_CD_ileum$qvalue<q_threshold_multi_trans),]
  tPWs_CD_rectum_sign <- tPWs_CD_rectum[(tPWs_CD_rectum$p.adjust<p_threshold_multi_trans)&(tPWs_CD_rectum$qvalue<q_threshold_multi_trans),]
  tPWs_UC_ileum_sign <- tPWs_UC_ileum[(tPWs_UC_ileum$p.adjust<p_threshold_multi_trans)&(tPWs_UC_ileum$qvalue<q_threshold_multi_trans),]
  tPWs_UC_rectum_sign <- tPWs_UC_rectum[(tPWs_UC_rectum$p.adjust<p_threshold_multi_trans)&(tPWs_UC_rectum$qvalue<q_threshold_multi_trans),]
  
  mPWs_CD_sign <- mPWs_CD[(mPWs_CD$probabilities< p_threshold_multi_mets)&(mPWs_CD$HMDBsInPWs>nMetsPathway)&(mPWs_CD$ProteinsInPWs>nProteinsPathway),]
  mPWs_UC_sign <- mPWs_UC[(mPWs_UC$probabilities< p_threshold_multi_mets)&(mPWs_UC$HMDBsInPWs>nMetsPathway)&(mPWs_UC$ProteinsInPWs>nProteinsPathway),]
  
  # Determine overlap
  overlapCD_all <- intersect(intersect(tPWs_CD_ileum_sign[,1], tPWs_CD_rectum_sign[,1]), mPWs_CD_sign[,1])
  overlapCD_ileum <- intersect(tPWs_CD_ileum_sign[,1], mPWs_CD_sign[,1])
  overlapCD_rectum <- intersect(tPWs_CD_rectum_sign[,1], mPWs_CD_sign[,1])
  
  overlapUC_all <- intersect(intersect(tPWs_UC_ileum_sign[,1], tPWs_UC_rectum_sign[,1]), mPWs_UC_sign[,1])
  overlapUC_ileum <- intersect(tPWs_UC_ileum_sign[,1], mPWs_UC_sign[,1])
  overlapUC_rectum <- intersect(tPWs_UC_rectum_sign[,1], mPWs_UC_sign[,1])
  
  overlapBoth_all <- intersect(overlapCD_all,overlapUC_all)
  overlapBoth_ileum <- intersect(overlapCD_ileum, overlapUC_ileum)
  overlapBoth_rectum <- intersect(overlapCD_rectum, overlapUC_rectum)
  
  # Make final table
  finalTable <- unique(rbind.data.frame(mPWs_CD, mPWs_UC)[,c(1,2)])
  finalTable_CD_all <- finalTable[finalTable[,1] %in% overlapCD_all,]
  finalTable_CD_ileum <- finalTable[finalTable[,1] %in% overlapCD_ileum,]
  finalTable_CD_rectum <- finalTable[finalTable[,1] %in% overlapCD_rectum,]
  
  finalTable_UC_all <- finalTable[finalTable[,1] %in% overlapUC_all,]
  finalTable_UC_ileum <- finalTable[finalTable[,1] %in% overlapUC_ileum,]
  finalTable_UC_rectum <- finalTable[finalTable[,1] %in% overlapUC_rectum,]
  
  output <- list(
    finalTable_CD_all,
    finalTable_CD_ileum,
    finalTable_CD_rectum,
    finalTable_UC_all,
    finalTable_UC_ileum,
    finalTable_UC_rectum
  )
  names(output) <- c("CD_Both",
                     "CD_Ileum",
                     "CD_Rectum",
                     "UC_Both",
                     "UC_Ileum",
                     "UC_Rectum")
  
  return(output)
}

