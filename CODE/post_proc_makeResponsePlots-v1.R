
library(biomod2)
library(raster)
library(ggplot2)
library(dplyr)
library(ggExtra)

source("D:/MyDocs/R-dev/biomod2plus/R/biomod2plus.R")

load("./OUT/Allsprecords_ModObjects.RData")

setwd("./OUT")


DF <- values(current) %>% na.omit %>% as.data.frame
colnames(DF)

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##


allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
                 "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")

for(modAlgo in allModAlgos){
  cat("Processing algorithm:",modAlgo,".....\n\n")
  myMods <- BIOMOD_LoadModels(myBiomodModelOut, models=modAlgo)
  
  responsePlots(myBiomodModelOut, Data=DF, modelsToUse = myMods, 
                showVars="all", fixedVarMetric = 'mean', plotStdErr = TRUE,
                addMarginalPlot = FALSE,
                filePrefix=paste("ResponsePlot_",modAlgo,"_",sep=""), 
                outFolder = "./RespPlots", height=4, width=4, plot=FALSE, 
                save=TRUE)
  rm(list=myMods)
  cat(" done.\n\n")
}



rm(list=ls()[grepl("Allsprecords_",ls())])
