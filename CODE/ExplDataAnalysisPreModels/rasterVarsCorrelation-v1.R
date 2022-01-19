
library(raster)
library(dplyr)
library(tidyr)


selVarsClim <- c(1,11,12,18)
climVarsPath <- list.files("./DATA/RASTER/Clim", pattern=".tif$", full.names = TRUE)
bioVarNames <- paste("bio_",1:19,sep="")[selVarsClim]

forestVarsPath <- list.files("./DATA/RASTER/ForestVars", pattern=".tif$", full.names = TRUE)
rsVarsPath <- list.files("./DATA/RASTER/RemoteSensing", pattern=".tif$", full.names = TRUE)
soilVarsPath <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)
topoVarsPath <- list.files("./DATA/RASTER/Topo", pattern=".tif$", full.names = TRUE)

allVars <- c(forestVarsPath, rsVarsPath, soilVarsPath, topoVarsPath)
varNames <- c(bioVarNames, gsub(".tif","",basename(allVars)))

# rstStack <- raster::stack(raster::stack(climVarsPath)[[selVarsClim]],raster::stack(allVars))
# names(rstStack) <- varNames
# rstData <- na.omit(values(rstStack))
#
#saveRDS(rstData, file = "./DATA/fullRasterStack.rds")
rstData <- readRDS("./DATA/fullRasterStack.rds") %>% as.data.frame

corMatSpearman <- cor(rstData, method="spearman") %>% round(2)
corMatPearson <- cor(rstData, method="pearson") %>% round(2)

#write.csv(corMatSpearman,"./OUT/corMatSpearman.csv")
#write.csv(corMatPearson,"./OUT/corMatPearson.csv")

caret::findCorrelation(corMat, cutoff = 0.8, names=TRUE)
# Remove:
# [1] "EVPT_median_mb" "NDWI_median_mb" "clay_lucas_ib" 
# [4] "ndwi_IQR_mb"    "silt_lucas_ib"

caret::findCorrelation(corMat, cutoff = 0.7, names=TRUE)
# Remove: 
# [1] "EVPT_median_mb" "NDWI_median_mb" "evpt_IQR_mb"   
# [4] "soil_pH_ib"     "clay_lucas_ib"  "ndwi_IQR_mb"   
# [7] "silt_lucas_ib"