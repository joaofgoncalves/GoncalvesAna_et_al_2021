
library(ade4)
library(ape)
library(raster)
library(dplyr)
library(magrittr)


dirs <- list.dirs(recursive = FALSE,full.names = TRUE)[-1]
projList <- c("current")

spNames <- c("AllSpeciesRecords",
             "T_lusitanica",
             "T_ebejeri",
             "T_stenoptera",
             "T_cantabrica",
             "T_semiaptera",
             "T_iberica",
             "PRC")

pathLists_ens <- list(current=c())
pathLists_bin <- list(current=c())

for(i in 1:length(dirs)){
  
  spFolderName <- basename(dirs[i])
  
  for(j in 1:length(projList)){
    
    proj <- projList[j]
    pathLists_ens[[proj]] <- c(pathLists_ens[[proj]],
                               paste(dirs[i],"/proj_",proj,"/GeoTIFF/proj_",proj,"_",spFolderName,"_ensemble.tif", sep=""))
    
    pathLists_bin[[proj]] <- c(pathLists_bin[[proj]],
                               paste(dirs[i],"/proj_",proj,"/GeoTIFF/proj_",proj,"_",spFolderName,"_ensemble_TSSbin.tif", sep=""))
    
  }
}

## ---------------------------------------------------------------------- ##

for(proj in projList){
  
  assign(paste("rstEns_",proj,sep=""),
         stack(pathLists_ens[[proj]], bands=1))
  
  assign(paste("rstBin_",proj,sep=""),
         stack(pathLists_bin[[proj]], bands=1))
}

## ---------------------------------------------------------------------- ##

bin_current_DF <- rstBin_current %>% values %>% na.omit %>% `colnames<-`(spNames)

ens_current_DF <- rstEns_current %>% values %>% na.omit %>% `colnames<-`(spNames)


## ---------------------------------------------------------------------- ##



## Habitat suitability clustering ------------------------------

## ALL -----------------------
#d_ens_current <- as.dist(sqrt(2*(1 - abs(cor(ens_current_DF, method="pearson")))))
d_ens_current <- as.dist((1 - cor(ens_current_DF, method="spearman")) / 2)
hc_ens_current <- hclust(d_ens_current, "average")

png("HC_HabSuit_Spearman_All-v1.png", res=300, width=1500, height=2000)
plot(hc_ens_current, main = "All | Hab. suit. Spearman | Average")
dev.off()

## Species only----------------
d_ens_current <- as.dist((1 - cor(ens_current_DF[,-c(1,8)], method="spearman")) / 2)
hc_ens_current <- hclust(d_ens_current, "average")

png("HC_HabSuit_Spearman_SpOnly-v1.png", res=300, width=1500, height=2000)
plot(hc_ens_current, main = "Species only | Hab. suit. Spearman | Average")
dev.off()


## Binary data clustering --------------------------------------

## ALL --------------------------
# Jaccard distance
d_bin_current_jac <- ade4::dist.binary(t(bin_current_DF), method=1)
hc_bin_current_jac <- hclust(d_bin_current_jac, "average")

png("HC_Bin_Jacc_All-v1.png", res=300, width=1500, height=2000)
plot(hc_bin_current_jac)
dev.off()

# Sorensen/Dice distance
d_bin_current_sor <- ade4::dist.binary(t(bin_current_DF), method=5)
hc_bin_current_sor <- hclust(d_bin_current_sor, "average")

png("HC_Bin_Sorensen_All-v1.png", res=300, width=1500, height=2000)
plot(hc_bin_current_sor)
dev.off()

## Species only ------------------
# Jaccard distance
d_bin_current_jac <- ade4::dist.binary(t(bin_current_DF[,-c(1,8)]), method=1)
hc_bin_current_jac <- hclust(d_bin_current_jac, "average")

png("HC_Bin_Jacc_SpOnly-v1.png", res=300, width=1500, height=2000)
plot(hc_bin_current_jac, main = "Species only | Bin. Sorensen | Average")
dev.off()

# Sorensen/Dice distance
d_bin_current_sor <- ade4::dist.binary(t(bin_current_DF[,-c(1,8)]), method=5)
hc_bin_current_sor <- hclust(d_bin_current_sor, "average")

png("HC_Bin_Sorensen_SpOnly-v1.png", res=300, width=1500, height=2000)
plot(hc_bin_current_sor, main = "Species only | Bin. Sorensen | Average")
dev.off()



