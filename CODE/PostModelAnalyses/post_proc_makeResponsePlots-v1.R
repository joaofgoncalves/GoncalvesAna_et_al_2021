
library(biomod2)
library(raster)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(tidyr)
library(biomod2plus)

#setwd("./OUT")

#source("D:/R-dev/biomod2plus/R/biomod2plus.R")

objs <- list.files(pattern=".RData$")
  
#   c(
# "Allsprecords_ModObjects.RData",
# "AriasellaLusitanica_ModObjects.RData",
# "AriasellaSp1_ModObjects.RData",
# "Cluster1_ModObjects.RData",
# "Cluster2_ModObjects.RData",
# "Cluster4_ModObjects.RData"
#           )

outFolders <- c(
"Allsprecords",
"AriasellaLusitanica",
"AriasellaSp1",
"AriasellaSp2",
"AriasellaSp3",
"Cluster1",
"Cluster2",
"Cluster4")
          

      
i <- 0
for(obj in objs){
  
  i <- i+1
  outFolderName <- outFolders[i]
  #source("D:/MyDocs/R-dev/biomod2plus/R/biomod2plus.R")
  
  outFn <- paste("./_ResponsePlots/",outFolderName,"/",sep="")
  
  if(!dir.exists(outFn)){
    dir.create(outFn)
  }
  
  load(obj)
  
  DF <- values(current) %>% na.omit %>% as.data.frame
  colnames(DF)
  
  modEvalData <- get_evaluations(myBiomodModelOut)
  
  evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
  
  # Top 3 modelling algorithms (by average TSS score) used 
  # to make response curves
  #
  # bestModAlgos <- apply(evalDF.TSS,1,mean) %>% 
  #   round(2) %>% 
  #   sort(decreasing = TRUE) %>% 
  #   names %>% 
  #   `[`(1:3)
  
  bestModAlgos<-"GBM"
  
  ## ------------------------------------------------------------------------ ##
  ##
  # allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
  #                  "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")
  # 
  # allModAlgos <- c("GBM", "GAM", "RF")
  
  # sort(unique(unlist(lapply(strsplit(modelsToUse,"_"), 
  # FUN = function(x) `[`(x,length(x))))))
  
  
  for(modAlgo in bestModAlgos){
    cat("Processing algorithm:",modAlgo,".....\n\n")
    
    modelsToUse <- BIOMOD_LoadModels(myBiomodModelOut, models=modAlgo)
    
    responsePlots(biomodModelOut = myBiomodModelOut, 
                  Data           = DF, 
                  modelsToUse    = modelsToUse, 
                  showVars       = "all", 
                  fixedVarMetric = 'mean', 
                  plotStdErr     = TRUE,
                  addMarginalPlot = TRUE, 
                  marginPlotType  = "histogram",
                  filePrefix      = "ResponsePlot_", 
                  fileSuffix      = paste("_",modAlgo,sep=""),
                  outFolder       = paste("./_ResponsePlots/",outFolderName,"/",sep=""), 
                  height          = 4, 
                  width           = 4, 
                  plot            = FALSE, 
                  save            = TRUE,
                  saveCSV         = TRUE)
    
    rm(list=modelsToUse)
    cat(" done.\n\n")
  }

  #rm(list=ls())

}



#rm(list=ls()[grepl("Allsprecords_",ls())])


## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

evalMetricPlot(myBiomodModelOut, evalMetric = "TSS", save = TRUE)

varImportancePlot(myBiomodModelOut,by = "all", save = TRUE)



