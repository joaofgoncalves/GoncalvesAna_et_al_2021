
library(raster)
library(rgdal)

ib <- raster("./DATA/RASTER/Mask/mask_ib1.tif")

rst <- stack(list.files("./DATA/RASTER/Clim/WClim_Present_1km/GeoTIFF",full.names = TRUE, pattern=".tif$"))

rstIbCrop <- crop(rst, ib)
rstIb <- mask(rstIbCrop, ib)

`dataType<-`(rstIb,"FLT4S")

writeRaster(rstIb,"./DATA/RASTER/Clim/WClim_Present_1km/Iberia/wclim_bio_present.tif")
