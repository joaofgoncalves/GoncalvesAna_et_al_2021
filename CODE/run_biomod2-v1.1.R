


library(raster)
library(biomod2)
library(dplyr)
library(readxl)
library(sp)
library(rgdal)
library(magrittr)
library(tools)
library(stringr)



## -------------------------------------------------------------------------------------- ##
## Aux functions ----
## -------------------------------------------------------------------------------------- ##

convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
  for(f in fl){
    cat("Converting file:\n",f," ...\n",sep="")
    fn <- file_path_sans_ext(basename(f))
    raster::stack(f) %>%
      raster::writeRaster(filename = paste(outputFolder,"/",fn,".tif",sep=""))
    cat("done!\n\n")
  }
}

countDistinct <- function(x) length(unique(x))

'%!in%' <- function(x,y)!('%in%'(x,y))

abbrevNames <- function(x) paste(str_to_title(unlist(strsplit(gsub("\\.","",x),"\\ "))),collapse="")

## function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

## -------------------------------------------------------------------------------------- ##
## Input parameters ----
## -------------------------------------------------------------------------------------- ##

# Selected variable names by group
# CHANGE THIS IF NEEDED!!...
#
selVarsClim <- paste("bio_",c(1,11,12,18),sep="")
selVarsForest <- c("broadlf_cover", "conif_cover", "tree_density")
selVarsRemoteSensing <- c("evi_IQR_mb","EVI_median_mb")
selVarsSoil <- c("awc_lucas_ib", "bulk_dens_lucas_ib", "coarse_frag_lucas_ib",
                 "sand_lucas_ib", "STU_EU_DEPTH_ROOTS_ib", "STU_EU_T_OC_ib")
selVarsTopo <- c("twi_ib")
  
# Column used to filter presence points
colToFilterBy <- "clusters" # Or "species"

# Specis or cluster name to model
# If the option "AllSpRecords" is used then no filter is applied 
# and all records are considered for modelling

#spNameSelected <- "cluster 2"
spNameSelected <- "AllSpRecords"

# Names for the 'projective' raster stacks
# See the section below where the raster data is loaded
# Put the name of each raster stack intended for projection in here
projNames <- c("current")

lonlatColNames <- c("longitude","latitude")


## -------------------------------------------------------------------------------------- ##
## Read input data ----
## -------------------------------------------------------------------------------------- ##


#spRecords <- read_excel("C:\\Users\\AnaGoncalves\\Dropbox\\MSc_AGoncalves//DATA//TABLES//matriz_locais.xlsx") %>% as.data.frame
spRecords <- read_excel("./DATA/TABLES/matriz_locais.xlsx") %>% as.data.frame

# Filter records by species or cluster
if(spNameSelected!="AllSpRecords"){
  spRecordsSelected <- spRecords %>% filter(!!as.name(colToFilterBy)==spNameSelected)
  
  spRecordPoints <- SpatialPoints(spRecordsSelected[,lonlatColNames],
                                  proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
}else{ # Don't filter records --- use them all!! ;-)

  spRecordPoints <- SpatialPoints(spRecords[,lonlatColNames],
                                  proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
}


## -------------------------------------------------------------------------------------- ##
## Load environmental data by date/scenario ----
## -------------------------------------------------------------------------------------- ##

# List of raster files from which to load the environmental variables
# USE RELATIVE PATHS.. NOT ABSOLUTE ONES...
climVars <- "./DATA/RASTER/Clim/wclim_bio_present.tif"
forestVars <- list.files("./DATA/RASTER/ForestVars", pattern=".tif$", full.names = TRUE)
remoteSensingVars <- list.files("./DATA/RASTER/RemoteSensing", pattern=".tif$", full.names = TRUE) 
soilVars <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE) 
topoVars <- list.files("./DATA/RASTER/Topo", pattern=".tif$", full.names = TRUE) 

# List of names for each group of environmental variables
climVarNames <- paste("bio_",1:19,sep="")
forestVarNames <- file_path_sans_ext(basename(forestVars))
remoteSensingVarNames <- file_path_sans_ext(basename(remoteSensingVars))
soilVarNames <- file_path_sans_ext(basename(soilVars))
topoVarNames <- file_path_sans_ext(basename(topoVars))

#current <- stack("C:\\Users\\AnaGoncalves\\Dropbox\\MSc_AGoncalves//DATA//RASTER//Clim//wclim_bio_present.tif")
current <- stack(c(climVars, forestVars, remoteSensingVars, soilVars, topoVars))

# Define raster layer names
names(current) <- c(climVarNames, forestVarNames, remoteSensingVarNames, 
                    soilVarNames, topoVarNames)

# Agreggate the selected var names
selVarsAll <- c(selVarsClim, selVarsForest, selVarsRemoteSensing, 
                selVarsSoil, selVarsTopo)

# Make the final raster stack with all the selected variables
current <- current[[selVarsAll]]

current <- stack(mask(current, intersect_mask(current)))

# proj2050 <- stack("")
# names(proj2050) <- paste("bio_",1:19,sep="")
# proj2050 <- proj2050[[selVars]]


## -------------------------------------------------------------------------------------- ##
## Run biomod2 ----
## -------------------------------------------------------------------------------------- ##

#setwd("C:\\Users\\AnaGoncalves\\Dropbox\\MSc_AGoncalves//OUT")
setwd("./OUT")

sp <- abbrevNames(spNameSelected)

# Number of training points
Npresences <- length(spRecordPoints)

# Set up the biomod data object for calibration
# Number of PA sets = 10
# Number of PA's per set = Number of presences
# PA selection strategy = random
#
myBiomodData <- BIOMOD_FormatingData(resp.var = spRecordPoints,
                                     expl.var = current,
                                     resp.name = sp,
                                     PA.nb.rep = 3, # Number of pseudo-absences sets
                                     PA.nb.absences = ifelse(Npresences < 500, Npresences*10, Npresences), # Nr of pseudo-absences
                                     PA.strategy = 'random') # PA generation method

# Model hyperparameters
# GAM: changes k=4 to avoid overly complex models
#
myBiomodOptions <- BIOMOD_ModelingOptions(GAM = list(k = 3),
                                          MAXENT.Phillips = list(threshold=FALSE,
                                                                 hinge=FALSE),
                                          GBM = list(n.trees = 1500))
#print(myBiomodOptions)


## -------------------------------------------------------------------------------------- ##
## Calibrate models ----
## -------------------------------------------------------------------------------------- ##


myBiomodModelOut <- BIOMOD_Modeling(
  data = myBiomodData, # Input data
  models = c('GLM','GBM','GAM','CTA','ANN',
             'FDA','MARS','RF','MAXENT.Phillips', 
             'MAXENT.Tsuruoka'), # Models to run
  models.options = myBiomodOptions,
  NbRunEval = 20, # Number of Evaluation runs
  DataSplit = 80, # Train percentage
  Prevalence = 0.5, # Prevalence between 0 and 1
  VarImport = 5, # Nr of rounds to evaluate variables
  models.eval.meth = c('TSS','ROC','KAPPA'), # Evaluation metrics
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = TRUE) # Model with all data?


# Get model evaluation values
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

# Print ROC scores
print(myBiomodModelEval["ROC","Testing.data",,,])
print(myBiomodModelEval["TSS","Testing.data",,,])

# Get boxplot stats
print(fivenum(as.numeric(myBiomodModelEval["ROC","Testing.data",,,])))
print(fivenum(as.numeric(myBiomodModelEval["TSS","Testing.data",,,])))

# Save evaluation metrics from the arrays
evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
evalDF.KAPPA <- as.data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])

write.csv(evalDF.ROC, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_ROC.csv",sep=""))
write.csv(evalDF.TSS, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_TSS.csv",sep=""))
write.csv(evalDF.KAPPA, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_KAPPA.csv",sep=""))


# Calculate variable importance across all PA sets, eval rouns and algorithms 
varImportance <- get_variables_importance(myBiomodModelOut)
varImportanceByVariableAVG <- apply(varImportance,1,mean)
varImportanceByVariableSTD <- apply(varImportance,1,sd)
vimpDF <- data.frame(cnames=names(varImportanceByVariableAVG),
                     vimpAVG = varImportanceByVariableAVG, 
                     varImpSTD=varImportanceByVariableSTD) %>% 
          arrange(desc(vimpAVG))

write.csv(vimpDF, file = paste(getwd(),"/",sp,"/",sp,"_varImportance.csv",sep=""))


## -------------------------------------------------------------------------------------- ##
## Perform ensemble modelling ----
## -------------------------------------------------------------------------------------- ##

# Ensemble all partial models
# Use the top 5% models
# Use ROC as the threshold metric
# Ensemble methods are: mean, median and weighted mean

quantileThresh <- quantile(evalDF.TSS, probs=0.75, na.rm=TRUE)
if(quantileThresh==1) quantileThresh <- quantile(evalDF.TSS, probs=0.70, na.rm=TRUE)
if(quantileThresh==1) quantileThresh <- quantile(evalDF.TSS, probs=0.65, na.rm=TRUE)

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by = 'all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = quantileThresh,
                                      prob.mean = TRUE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      prob.median = FALSE,
                                      committee.averaging = FALSE,
                                      prob.mean.weight = FALSE,
                                      prob.mean.weight.decay = 'proportional')


# Get evaluation scores for the Ensemble Modelling stage
emEvalDF <- as.data.frame(get_evaluations(myBiomodEM))
write.csv(emEvalDF, file = paste(getwd(),"/",sp,"/",sp,"_EnsMod_evalDF_AllMetrics.csv",sep=""))


## -------------------------------------------------------------------------------------- ##
## Obtain spatiotemporal projections ----
## -------------------------------------------------------------------------------------- ##

# Models to consider in the ensemble and projection
modelsToUse <- get_kept_models(myBiomodEM, 1)



for(projName in projNames){
  
  # Obtain spatiotemporal projections
  myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env = get(projName),
                                    proj.name = projName, ## Name of the projection from above variable proj.name
                                    selected.models = modelsToUse,
                                    filtered.meth = NULL,
                                    binary.meth = NULL,
                                    compress = 'gzip',
                                    clamping.mask = TRUE,
                                    output.format = '.grd',
                                    do.stack = TRUE)
  
  
  # Perform the ensembling of projections
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           binary.meth = c('TSS','ROC','KAPPA'),
                                           EM.output = myBiomodEM,
                                           output.format = '.grd')
  
  # Convert all output raster files to GeoTIFF
  inFolder <- paste(getwd(),"/",sp,"/proj_",projName,sep="")
  outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
  dir.create(outFolder)
  
  convertToGeoTIFF(inFolder, outFolder)
  
} 

save.image(file=paste(sp,"ModObjects.RData",sep="_"))

