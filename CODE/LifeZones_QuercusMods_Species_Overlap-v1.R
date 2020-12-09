
library(raster)
library(dplyr)
library(magrittr)

setwd("D:/Projects/Colab_AGoncalves/MSc_AGoncalves/MODS")


lifeZones <- raster("D:/R-dev/QuercusSDMs/OUT_/LifeZones/Biogeo_regions_life_zones_int_v1.tif")

lifeZoneClasses <- c(
  "Mediterranean",
  "Temperate_Eurosiberian",
  "Submediterranean",
  "Eurosiberian_SubmediterraneanTransition"
)

dirs <- list.dirs(recursive = FALSE,full.names = TRUE)[-c(5,6)]
projList <- c("current")

spNames <- c("AllSpeciesRecords",
             "T_lusitanica",
             "T_ebejeri",
             #"T_stenoptera",
             #"T_cantabrica",
             "T_semiaptera",
             "T_iberica",
             "PRC")

pathLists_bin <- list(current=c())
pathLists_ens <- list(current=c())

for(i in 1:length(dirs)){
  
  spFolderName <- basename(dirs[i])
  
  for(j in 1:length(projList)){
    
    proj <- projList[j]
    
    pathLists_ens[[proj]] <- c(pathLists_ens[[proj]],
                               paste(dirs[i],"/proj_",proj,"/GeoTIFF/proj_",proj,"_",
                                     spFolderName,"_ensemble.tif", sep=""))
    
    
    pathLists_bin[[proj]] <- c(pathLists_bin[[proj]],
                               paste(dirs[i],"/proj_",proj,"/GeoTIFF/proj_",proj,"_",
                                     spFolderName,"_ensemble_TSSbin.tif", sep=""))
    
  }
}

## ---------------------------------------------------------------------- ##

for(proj in projList){
  
  assign(paste("rstEns_",proj,sep=""),
         stack(pathLists_ens[[proj]], bands=1))
  
  assign(paste("rstBin_",proj,sep=""),
         stack(pathLists_bin[[proj]], bands=1))
}

print(rstBin_current)
print(rstEns_current)

names(rstBin_current) <- spNames
names(rstEns_current) <- spNames


## ---------------------------------------------------------------------- ##


fl1 <- list.files("D:/R-dev/QuercusSDMs/OUT_/FINAL_MODS_1Km_CVV_et_al_2020_IJGI/CLIM_EFA/BIN_MAPS", 
                  pattern=".tif$",full.names = TRUE)

qspNames <- substr(basename(fl1), 1, 4)

quercusMods <- stack(fl1)
names(quercusMods) <- qspNames


fl2 <- list.files("D:/R-dev/QuercusSDMs/OUT_/FINAL_MODS_1Km_CVV_et_al_2020_IJGI/CLIM_EFA/HAB_SUIT", 
                  pattern=".tif$",full.names = TRUE)

qspNames <- substr(basename(fl2), 1, 4)

quercusModsHS <- stack(fl2)
names(quercusModsHS) <- qspNames

## ---------------------------------------------------------------------- ##


iouMatrixLZ <- matrix(NA,nrow=length(spNames), ncol=length(lifeZoneClasses))
colnames(iouMatrixLZ) <- lifeZoneClasses
rownames(iouMatrixLZ) <- spNames

intCovMatrixLZ <- matrix(NA,nrow=length(spNames), ncol=length(lifeZoneClasses))
colnames(intCovMatrixLZ) <- lifeZoneClasses
rownames(intCovMatrixLZ) <- spNames

iouMatrix <- matrix(NA,nrow=length(spNames), ncol=length(qspNames))
colnames(iouMatrix) <- qspNames
rownames(iouMatrix) <- spNames

intCovMatrix <- matrix(NA,nrow=length(spNames), ncol=length(qspNames))
colnames(intCovMatrix) <- qspNames
rownames(intCovMatrix) <- spNames

corMatrix <- matrix(NA,nrow=length(spNames), ncol=length(qspNames))
colnames(corMatrix) <- qspNames
rownames(corMatrix) <- spNames
  
pb <- txtProgressBar(min=1,max=length(spNames)*length(qspNames),style=3)




i <- 0
for(spName in spNames){
  
  sp_s <- sum(na.omit(as.numeric(values(rstBin_current[[spName]]))))
  
  for(cl in c(1:4)){

      int <- rstBin_current[[spName]] & (lifeZones==cl)
      uni <- rstBin_current[[spName]] | (lifeZones==cl)
      
      int_s <- sum(na.omit(as.numeric(values(int))))
      uni_s <- sum(na.omit(as.numeric(values(uni))))
      
      
      iouMatrixLZ[spName, cl] <- int_s / uni_s
      intCovMatrixLZ[spName, cl] <- int_s / sp_s
  }
  
  
  for(qspName in qspNames){
    i<-i+1
    int <- rstBin_current[[spName]] & quercusMods[[qspName]]
    uni <- rstBin_current[[spName]] | quercusMods[[qspName]]
    
    int_s <- sum(na.omit(as.numeric(values(int))))
    uni_s <- sum(na.omit(as.numeric(values(uni))))
    
    iouMatrix[spName, qspName] <- int_s / uni_s
    intCovMatrix[spName, qspName] <- int_s / sp_s
    
    corValue <- stack(rstEns_current[[spName]], quercusModsHS[[qspName]]) %>% 
      values %>% na.omit %>% cor(method="spearman") %>% `[`(2,1)
    
    corMatrix[spName, qspName] <- corValue
    
    setTxtProgressBar(pb, i)
    
  }
  
}

iouMatrix %>% round(2)
corMatrix %>% round(2)
iouMatrixLZ %>% round(2)

intCovMatrix %>% round(2)
intCovMatrixLZ %>% round(2)
