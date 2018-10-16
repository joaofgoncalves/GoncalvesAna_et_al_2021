


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
## Don't change this part!!... ;-)
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


## -------------------------------------------------------------------------------------- ##
## Input parameters ----
## -------------------------------------------------------------------------------------- ##

# Selected variables for bioclimatic indices
selVarsClim <- c(1,11,12,18)

# Column used to filter presence points
colToFilterBy <- "clusters" # Or "clusters"

# Specis or cluster name to model
spNameSelected <- "cluster 4"

# Names for the 'projective' raster stacks
# See the section below where the raster data is loaded
# Put the name of each raster stack intended for projection in here
projNames <- c("current")

# Path to the input table containing species presence records
inputSpRecordsTablePath <- "./DATA/TABLES/matriz_locais.xlsx"

# Column names for longitude and latitude
lonlatColNames <- c("longitude","latitude")

# Output folder location (all biomod2 outputs will be placed here)
outputFolder <- "./OUT"

## -------------------------------------------------------------------------------------- ##
## Read input data ----
## -------------------------------------------------------------------------------------- ##

# Load the input table with species presences
spRecords <- read_excel(inputSpRecordsTablePath) %>% as.data.frame

# Subset the table to include only the target species
spRecordsSelected <- spRecords %>% filter(!!as.name(colToFilterBy)==spNameSelected)

# Make a geographical objected with the selected points
spRecordPoints <- SpatialPoints(spRecordsSelected[,lonlatColNames],
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))


## -------------------------------------------------------------------------------------- ##
## Load environmental data by date/scenario ----
## -------------------------------------------------------------------------------------- ##

# Load predictive variables in raster format using raster::stack
#
current <- stack("./DATA/RASTER/Clim/wclim_bio_present.tif")
names(current) <- paste("bio_",1:19,sep="")
current <- current[[selVarsClim]]

# proj2050 <- stack("")
# names(proj2050) <- paste("bio_",1:19,sep="")
# proj2050 <- proj2050[[selVars]]


## -------------------------------------------------------------------------------------- ##
## Run biomod2 ----
## -------------------------------------------------------------------------------------- ##

# Set the working directory where outputs will be placed
setwd(outputFolder)

# Abbreviated species name
sp <- abbrevNames(spNameSelected)

# Number of training points
Npresences <- length(spRecordPoints)

# Set up the biomod data object for calibration
# Number of PA sets = 10
# Number of PA's per set = Number of presences
# PA selection strategy = random
#
myBiomodData <- BIOMOD_FormatingData(resp.var = spRecordPoints, # Input presence records
                                     expl.var = current, # Predictive variables as a raster stack
                                     resp.name = sp, # Name of the target-species
                                     PA.nb.rep = 1, # Number of pseudo-absences sets
                                     PA.nb.absences = ifelse(Npresences < 500, Npresences*10, Npresences), # Nr of pseudo-absences
                                     PA.strategy = 'random') # PA generation method

# Model hyperparameters
# GAM: changes k=4 to avoid overly complex models
#
myBiomodOptions <- BIOMOD_ModelingOptions(GAM = list(k = 4),
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
  NbRunEval = 3, # Number of Evaluation run
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
vimpDF <- data.frame(vimpAVG = varImportanceByVariableAVG, varImpSTD=varImportanceByVariableSTD)

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

