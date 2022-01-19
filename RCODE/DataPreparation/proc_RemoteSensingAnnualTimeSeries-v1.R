
library(raster)

mask <- raster("./DATA/RASTER/Mask/mask_ib1.tif")

varPaths <- list.files("C:/Users/JG/Downloads/rs_vars", pattern = ".tif$", full.names = TRUE)


for(i in 1:length(varPaths)){
  
  # Load the time series for each indicator
  rst <- stack(varPaths[i])
  
  # Calculate the multi-annual average between 2008-2017 (10yrs)
  rstMean <- calc(rst[[8:17]], fun = mean)
  
  # Re-project to 1km
  rstMeanProj <- projectRaster(rstMean, to = mask)
  
  # Save/write data
  writeRaster(rstMeanProj,paste("./DATA/RASTER/RemoteSensing/",basename(varPaths[i]),sep=""),
              overwrite=TRUE)
  
  cat("-> Finished file:",basename(varPaths[i]),"......\n\n")
}


