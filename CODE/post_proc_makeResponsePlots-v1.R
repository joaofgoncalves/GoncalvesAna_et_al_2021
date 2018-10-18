
library(biomod2)
library(raster)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(tidyr)
library(biomod2plus)
#source("D:/MyDocs/R-dev/biomod2plus/R/biomod2plus.R")

setwd("./OUT")

load("Allsprecords_ModObjects.RData")


DF <- values(current) %>% na.omit %>% as.data.frame
colnames(DF)


## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

# allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
#                  "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")
# 
allModAlgos <- c("GBM", "GAM", "RF")

for(modAlgo in allModAlgos){
  cat("Processing algorithm:",modAlgo,".....\n\n")
  myMods <- BIOMOD_LoadModels(myBiomodModelOut, models=modAlgo)
  
  responsePlots(myBiomodModelOut, Data=DF, modelsToUse = myMods, 
                showVars="all", fixedVarMetric = 'mean', plotStdErr = TRUE,
                addMarginalPlot = TRUE, marginPlotType = "histogram",
                filePrefix="ResponsePlot_", fileSuffix=paste("_",modAlgo,sep=""),
                outFolder = "./Allsprecords/ResponsePlots", height=4, 
                width=4, plot=FALSE, save=TRUE)
  rm(list=myMods)
  cat(" done.\n\n")
}



rm(list=ls()[grepl("Allsprecords_",ls())])


## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

evalMetricPlot(myBiomodModelOut, evalMetric = "TSS", save = TRUE)

varImportancePlot(myBiomodModelOut,by = "all", save = TRUE)



