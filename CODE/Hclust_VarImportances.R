
library(dplyr)
library(tidyr)

baseFolder <- "D:/Projects/Colab_AGoncalves/MSc_AGoncalves/MODS"

dirList <- list.dirs(baseFolder,recursive = FALSE,full.names = FALSE)[-1]

spNames <- c("AllSpeciesRecords",
             "T_lusitanica",
             "T_ebejeri",
             "T_stenoptera",
             "T_cantabrica",
             "T_semiaptera",
             "T_iberica",
             "PRC")
i <- 0
for(dirName in dirList){
  
  i <- i + 1 
  varImpCSVfile <- paste(baseFolder,"/",dirName,"/",dirName,"_varImportance.csv",sep="")
  tmp <- data.frame(spName = spNames[i], read.csv(varImpCSVfile)[,-c(1,4)])
  
  if(i==1){
    varImpBySp <- tmp
  }else{
    varImpBySp <- bind_rows(varImpBySp, tmp)
  }
  
}

varImpBySp2 <- pivot_wider(varImpBySp, 
            names_from = "spName",
            id_cols = c("cnames"), 
            values_from = "vimpAVG")


#d <- as.dist(1 - abs(cor(varImpBySp2[,-1], method = "spearman")))
d <- as.dist((1 - cor(varImpBySp2[,-1], method = "spearman")) / 2)

#d <- dist(varImpBySp2[,-1])

hc.wardD2 <- hclust(d, method="ward.D2")
hc.avg <- hclust(d, method="average")

png(filename="HClust_VarImp_All.png",width = 2800, height = 2000, res = 300)
par(mfrow=c(1,2))
plot(hc.wardD2, main = "HClust Var. imp. | Ward D2")
plot(hc.avg, main = "HClust Var. imp. | Average")
dev.off()

varImpBySp3 <- pivot_wider(varImpBySp %>% filter(!(spName %in% c("AllSpeciesRecords","PRC"))), 
                           names_from = "spName",
                           id_cols = c("cnames"), 
                           values_from = "vimpAVG")


#d <- as.dist(1 - abs(cor(varImpBySp3[,-1], method = "spearman")))
d <- as.dist((1 - cor(varImpBySp3[,-1], method = "spearman")) / 2)
#d <- dist(varImpBySp3[,-1])

hc.wardD2 <- hclust(d, method="ward.D2")
hc.avg <- hclust(d, method="average")

png(filename="HClust_VarImp_SpeciesOnly.png",width = 2800, height = 2000, res = 300)
par(mfrow=c(1,2))
plot(hc.wardD2, main = "HClust Var. imp. | Ward D2")
plot(hc.avg, main = "HClust Var. imp. | Average")
dev.off()




