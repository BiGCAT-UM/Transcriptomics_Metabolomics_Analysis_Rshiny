
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
htxCount <- read.csv(paste0(WORK.DIR,"/1-data_preprocessing/htxCount.csv"))
sampleLabels <- read.csv(paste0(WORK.DIR,"/1-data_preprocessing/sampleLabels.csv"),row.names = 1 )

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
mappingTranscriptomics <- function (){
  
  #input data is all files starts with "table" under the folder 2-differential_gene_expression_analysis/statsmodel/
  #output files will be in "3-identifier_mapping" folder with the same structure IDMapping_CD OR IDMapping_UC
  
}

pathwayAnalysisTranscriptomics <- function(){
  
}

createHeatmap <- function(){
  
}

networkAnalysis <- function(){
  
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
  
  # png(paste("FC_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
  # hist(toptab[,"Fold Change"],main=paste("adapted fold change histogram for",colnames(design)[i]),
  #      xlab="adapted fold changes",col="green3",breaks=120,cex.axis=1.2,cex.lab=1.2)
  # dev.off()
  
  ## Visualize the data after the transformation (for one sample to get an idea of suitability of transformation:
  #create histogram for original distribution for first column with data
  png(paste("7-metabolite_data_preprocessing/normalized/CD_histogram_raw",".png"),width=1000,height=1000)
  hist(mSet_CD[,3], col='steelblue', main='Original')
  dev.off()
  cat("raw CD histogram created")
  
  #create histogram for log-transformed distribution 
  png(paste("7-metabolite_data_preprocessing/normalized/CD_histogram_norm",".png"),width=1000,height=1000)
  hist(mSet_transformed[,3], col='steelblue', main='Original')
  dev.off()
  cat("norm CD histogram created")
  
  
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
  png(paste("7-metabolite_data_preprocessing/normalized/UC_histogram_raw",".png"),width=1000,height=1000)
  hist(mSet_UC[,3], col='steelblue', main='Original')
  dev.off()
  cat("raw UC histogram created")
  
  #create histogram for log-transformed distribution 
  png(paste("7-metabolite_data_preprocessing/normalized/UC_histogram_norm",".png"),width=1000,height=1000)
  hist(mSet_transformed[,3], col='steelblue', main='Original')
  dev.off()
  cat("norm UC histogram created")
  
  
  
}

####################################################STATISTICAL ANALYSIS FUNCTIONS  #########################################

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