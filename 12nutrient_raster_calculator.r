####################################################################################
# Calculate Soil Nutrient Availability globally for static and seasonal effects
# Author: Jared Streich
# Version 1.1.0
# Email: ju0@ornl.gov, if not at ORNL streich.jared@gmail.com
# Date: 2021-06-15 and 16th
####################################################################################


####################################################################################
################################ Load Libraries ####################################
####################################################################################
library(raster)
library(dplyr)


####################################################################################
################################ Start Script ######################################
####################################################################################

##### Set color palette
cols <- colorRampPalette(c("grey30", "cadetblue3", "deepskyblue1", "olivedrab2", "sienna3", "red3"))( 200 )


##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- ph
bkgrnd[bkgrnd < 0] <- 0
bkgrnd[bkgrnd >= 0] <- 0
writeRaster(bkgrnd, filename = "bkgrnd_2021-06-15.tif", format = "GTiff", overwrite = T)
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
bkgrnd
ph

##### Nitrogen Availability
nph <- ph
nph[nph >= 6 & nph < 6.5] <- 6.0
nph[nph >= 0 & nph < 4] <- 0.01
nph[nph >= 4 & nph < 4.5] <- 1.2
nph[nph >= 4.5 & nph < 5] <- 3.0
nph[nph >= 5 & nph < 5.5] <- 4.0
nph[nph >= 5.5 & nph < 6] <- 4.5
nph[nph >= 6.5 & nph < 7] <- 6.0
nph[nph >= 7 & nph < 7.5] <- 6.0
nph[nph >= 7.5 & nph < 8] <- 6.0
nph[nph >= 8 & nph < 8.5] <- 4.0
nph[nph >= 8.5 & nph < 9] <- 3.5
nph[nph >= 9 & nph < 9.5] <- 2.6
nph[nph >= 9.5 & nph < 10] <- 1.3
nph[nph >= 10 & nph < 20] <- 0.01
nph <- nph/6
nph[nph >= 1] <- 0.5
nph.0 <- focal(nph, w = matrix(1,11,11), mean, na.rm = T)
nph.0 <- nph.0 + bkgrnd
writeRaster(nph.0, filename = "Nit_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(nph, filename = "Nit_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)
plot(nph, col = cols)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
"nitrogen layer made"

##### Phosphorus Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
bkgrnd
ph

pph <- ph
pph[pph >= 6 & pph < 6.5] <- 6
pph[pph >= 0 & pph < 4] <- 0.01
pph[pph >= 4 & pph < 4.5] <- 0.1
pph[pph >= 4.5 & pph < 5] <- 0.67
pph[pph >= 5 & pph < 5.5] <- 1.98
pph[pph >= 5.5 & pph < 6] <- 4.2
pph[pph >= 6.5 & pph < 7] <- 6
pph[pph >= 7 & pph < 7.5] <- 6
pph[pph >= 7.5 & pph < 8] <- 6.0
pph[pph >= 8 & pph < 8.5] <- 4.1
pph[pph >= 8.5 & pph < 9] <- 2.24
pph[pph >= 9 & pph < 9.5] <- 6
pph[pph >= 9.5 & pph < 10] <- 6
pph[pph >= 10 & pph < 20] <- 6
pph <- pph/6
pph[pph >= 1] <- 0.5
plot(pph, col = cols)
pph.0 <- focal(pph, w = matrix(1,11,11), mean, na.rm = T)
pph.0 <- pph.0 + bkgrnd
writeRaster(pph.0, filename = "phosph_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(pph, filename = "phosph_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
"phosphorus layer made"


##### Potasium Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

potph <- ph
potph[potph >= 6.5 & potph < 7] <- 6
potph[potph >= 0 & potph < 4] <- 0.01
potph[potph >= 4 & potph < 4.5] <- 0.1
potph[potph >= 4.5 & potph < 5] <- 1.54
potph[potph >= 5 & potph < 5.5] <- 3.2
potph[potph >= 5.5 & potph < 6] <- 4.6
potph[potph >= 6 & potph < 6.5] <- 6
potph[potph >= 7 & potph < 7.5] <- 6
potph[potph >= 7.5 & potph < 8] <- 6.0
potph[potph >= 8 & potph < 8.5] <- 6
potph[potph >= 8.5 & potph < 9] <- 6
potph[potph >= 9 & potph < 9.5] <- 6
potph[potph >= 9.5 & potph < 10] <- 6
potph[potph >= 10 & potph < 20] <- 6
potph <- potph/6
potph[potph >= 1] <- 0.5
plot(potph, col = cols)
potph.0 <- focal(potph, w = matrix(1,11,11), mean, na.rm = T)
potph.0 <- potph.0 + bkgrnd
writeRaster(potph.0, filename = "Potasium_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(potph, filename = "Potasium_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
"potasium layer made"

##### Sulphur Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

sph <- ph
sph[sph >= 6 & sph < 6.5] <- 5.7
sph[sph >= 0 & sph < 4] <- 0.01
sph[sph >= 4 & sph < 4.5] <- 0.1
sph[sph >= 4.5 & sph < 5] <- 1.54
sph[sph >= 5 & sph < 5.5] <- 3.2
sph[sph >= 5.5 & sph < 6] <- 4.3
sph[sph >= 6.5 & sph < 7] <- 6
sph[sph >= 7 & sph < 7.5] <- 6
sph[sph >= 7.5 & sph < 8] <- 6.0
sph[sph >= 8 & sph < 8.5] <- 6
sph[sph >= 8.5 & sph < 9] <- 6
sph[sph >= 9 & sph < 9.5] <- 6
sph[sph >= 9.5 & sph < 10] <- 6
sph[sph >= 10 & sph < 20] <- 6
sph <- sph/6
sph[sph >= 1] <- 0.5
sph.0 <- focal(sph, w = matrix(1,11,11), mean, na.rm = T)
sph.0 <- sph.0 + bkgrnd
writeRaster(sph.0, filename = "Sulph_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(sph, filename = "Sulph_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
"Sulphur layer made"

##### Calcium Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

calph <- ph
calph[calph >= 6 & calph < 6.5] <- 3.8
calph[calph >= 0 & calph < 4] <- 0.01
calph[calph >= 4 & calph < 4.5] <- 0.1
calph[calph >= 4.5 & calph < 5] <- 1.1
calph[calph >= 5 & calph < 5.5] <- 1.93
calph[calph >= 5.5 & calph < 6] <- 2.7
calph[calph >= 6.5 & calph < 7] <- 4.6
calph[calph >= 7 & calph < 7.5] <- 6
calph[calph >= 7.5 & calph < 8] <- 6
calph[calph >= 8 & calph < 8.5] <- 5.5
calph[calph >= 8.5 & calph < 9] <- 4.9
calph[calph >= 9 & calph < 9.5] <- 3.0
calph[calph >= 9.5 & calph < 10] <- 2.6
calph[calph >= 10 & calph < 20] <- 1.0
calph <- calph/6
calph[calph >= 1] <- 0.5
calph.0 <- focal(calph, w = matrix(1,11,11), mean, na.rm = T)
calph.0 <- calph.0 + bkgrnd
writeRaster(calph.0, filename = "Cal_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(calph, filename = "Cal_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)
##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
"Calcium layer made"



##### Magnesium Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

mgph <- ph
mgph[mgph >= 6 & mgph < 6.5] <- 3.2
mgph[mgph >= 0 & mgph < 4] <- 0.01
mgph[mgph >= 4 & mgph < 4.5] <- 0.1
mgph[mgph >= 4.5 & mgph < 5] <- 0.5
mgph[mgph >= 5 & mgph < 5.5] <- 1.73
mgph[mgph >= 5.5 & mgph < 6] <- 2.5
mgph[mgph >= 6.5 & mgph < 7] <- 4.0
mgph[mgph >= 7 & mgph < 7.5] <- 5
mgph[mgph >= 7.5 & mgph < 8] <- 5.6
mgph[mgph >= 8 & mgph < 8.5] <- 5.5
mgph[mgph >= 8.5 & mgph < 9] <- 5.4
mgph[mgph >= 9 & mgph < 9.5] <- 4.67
mgph[mgph >= 9.5 & mgph < 10] <- 2.8
mgph[mgph >= 10 & mgph < 20] <- 1.35
mgph <- mgph/5.6
mgph[mgph >= 1] <- 0.5
mgph.0 <- focal(mgph, w = matrix(1,11,11), mean, na.rm = T)
mgph.0 <- mgph.0 + bkgrnd
writeRaster(mgph.0, filename = "Magnesium_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(mgph, filename = "Magnesium_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Magnesium layer made"



##### Iron Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

iph <- ph
iph[iph >= 6 & iph < 6.5] <- 6
iph[iph >= 0 & iph < 4] <- 6
iph[iph >= 4 & iph < 4.5] <- 6
iph[iph >= 4.5 & iph < 5] <- 6
iph[iph >= 5 & iph < 5.5] <- 6
iph[iph >= 5.5 & iph < 6] <- 6
iph[iph >= 6.5 & iph < 7] <- 4.8
iph[iph >= 7 & iph < 7.5] <- 4.6
iph[iph >= 7.5 & iph < 8] <- 3.75
iph[iph >= 8 & iph < 8.5] <- 2.6
iph[iph >= 8.5 & iph < 9] <- 2
iph[iph >= 9 & iph < 9.5] <- 2
iph[iph >= 9.5 & iph < 10] <- 2
iph[iph >= 10 & iph < 20] <- 2
iph <- iph/6
iph[iph >= 1] <- 0.5
iph.0 <- focal(iph, w = matrix(1,11,11), mean, na.rm = T)
iph.0 <- iph.0 + bkgrnd
writeRaster(iph.0, filename = "Iron_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(iph, filename = "Iron_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Iron layer made"




##### Manganese Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

mnph <- ph
mnph[mnph >= 6 & mnph < 6.5] <- 6
mnph[mnph >= 0 & mnph < 4] <- 1
mnph[mnph >= 4 & mnph < 4.5] <- 1.7
mnph[mnph >= 4.5 & mnph < 5] <- 3.4
mnph[mnph >= 5 & mnph < 5.5] <- 6
mnph[mnph >= 5.5 & mnph < 6] <- 6
mnph[mnph >= 6.5 & mnph < 7] <- 6
mnph[mnph >= 7 & mnph < 7.5] <- 5.1
mnph[mnph >= 7.5 & mnph < 8] <- 4.6
mnph[mnph >= 8 & mnph < 8.5] <- 3.92
mnph[mnph >= 8.5 & mnph < 9] <- 2.6
mnph[mnph >= 9 & mnph < 9.5] <- 2.3
mnph[mnph >= 9.5 & mnph < 10] <- 1.7
mnph[mnph >= 10 & mnph < 20] <- 0.8
mnph <- mnph/6
mnph[mnph >= 1] <- 0.5
mnph.0 <- focal(mnph, w = matrix(1,11,11), mean, na.rm = T)
mnph.0 <- mnph.0 + bkgrnd
writeRaster(mnph.0, filename = "Manganese_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(mnph, filename = "Manganese_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Manganese layer made"



##### Boron Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

bph <- ph
bph[bph >= 6 & bph < 6.5] <- 6
bph[bph >= 0 & bph < 4] <- 1
bph[bph >= 4 & bph < 4.5] <- 2.2
bph[bph >= 4.5 & bph < 5] <- 3.4
bph[bph >= 5 & bph < 5.5] <- 3.7
bph[bph >= 5.5 & bph < 6] <- 5.25
bph[bph >= 6.5 & bph < 7] <- 6
bph[bph >= 7 & bph < 7.5] <- 6
bph[bph >= 7.5 & bph < 8] <- 4.1
bph[bph >= 8 & bph < 8.5] <- 2.8
bph[bph >= 8.5 & bph < 9] <- 1.5
bph[bph >= 9 & bph < 9.5] <- 2.3
bph[bph >= 9.5 & bph < 10] <- 5.5
bph[bph >= 10 & bph < 20] <- 5.5
bph <- bph/6
bph[bph >= 1] <- 0.5
bph.0 <- focal(bph, w = matrix(1,11,11), mean, na.rm = T)
bph.0 <- bph.0 + bkgrnd
writeRaster(bph.0, filename = "Boron_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(bph, filename = "Boron_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Boron layer made"



##### Molybdenum Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

molph <- ph
molph[molph >= 6 & molph < 6.5] <- 4
molph[molph >= 0 & molph < 4] <- 0.01
molph[molph >= 4 & molph < 4.5] <- 0.01
molph[molph >= 4.5 & molph < 5] <- 0.6
molph[molph >= 5 & molph < 5.5] <- 1.9
molph[molph >= 5.5 & molph < 6] <- 2.75
molph[molph >= 6.5 & molph < 7] <- 5.4
molph[molph >= 7 & molph < 7.5] <- 5.9
molph[molph >= 7.5 & molph < 8] <- 6
molph[molph >= 8 & molph < 8.5] <- 6
molph[molph >= 8.5 & molph < 9] <- 6
molph[molph >= 9 & molph < 9.5] <- 6
molph[molph >= 9.5 & molph < 10] <- 6
molph[molph >= 10 & molph < 20] <- 6
molph <- molph/6
molph[molph >= 1] <- 0.5
molph.0 <- focal(molph, w = matrix(1,11,11), mean, na.rm = T)
molph.0 <- molph.0 + bkgrnd
writeRaster(molph.0, filename = "Moleb_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(molph, filename = "Moleb_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Molybdenum layer made"



##### Copper/Zinc Availability
##### Read in pH data
ph <- raster("/Volumes/Smithers/Rasters/static_rasters/soil_pH_0-14.tif")
bkgrnd <- raster("bkgrnd_2021-06-15.tif")
ph

czph <- ph
czph[czph >= 6 & czph < 6.5] <- 6
czph[czph >= 0 & czph < 4] <- 0.01
czph[czph >= 4 & czph < 4.5] <- 1.2
czph[czph >= 4.5 & czph < 5] <- 2.9
czph[czph >= 5 & czph < 5.5] <- 6
czph[czph >= 5.5 & czph < 6] <- 6
czph[czph >= 6.5 & czph < 7] <- 6
czph[czph >= 7 & czph < 7.5] <- 6
czph[czph >= 7.5 & czph < 8] <- 4.4
czph[czph >= 8 & czph < 8.5] <- 3
czph[czph >= 8.5 & czph < 9] <- 1.5
czph[czph >= 9 & czph < 9.5] <- 1.5
czph[czph >= 9.5 & czph < 10] <- 1.5
czph[czph >= 10 & czph < 20] <- 1.5
czph <- czph/6
czph[czph >= 1] <- 0.5
czph.0 <- focal(czph, w = matrix(1,11,11), mean, na.rm = T)
czph.0 <- czph.0 + bkgrnd
writeRaster(czph.0, filename = "Copp-Zinc_soil_pH_2021-06-15_wIntrp.tif", format = "GTiff", overwrite = T)
writeRaster(czph, filename = "Copp-Zinc_soil_pH_2021-06-15.tif", format = "GTiff", overwrite = T)

##### Clear temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

"Copper/Zinc layer made"







#####################################################################################
############################ Get seasonal Data Prepared #############################
#####################################################################################

##### Set directory
setwd("/Volumes/Smithers/Rasters/")

soilwater <- list.files(pattern = "swc_fr_*")
temp <- list.files(pattern = "wc2.0_30s_tavg_*")

for(i in 1:length(soilwater)){
  rast <- raster(soilwater[i])
  val.new <- maxValue(rast)
  if(i == 1){
    val.pre <- val.new
  }
  else{
    val.pre <- c(val.pre, val.new)
  }
}
soil.max <- max(val.pre)
print(soil.max)


for(i in 1:length(temp)){
  rast <- raster(temp[i])
  val.new <- maxValue(rast)
  val.min.new <- minValue(rast)
  if(i == 1){
    val.pre <- val.new
    val.min.pre <- val.min.new
  }
  else{
    val.pre <- c(val.pre, val.new)
    val.min.pre <- c(val.min.pre, val.min.new)
  }
}
temp.range <- max(val.pre) - min(val.min.pre)
temp.range






i <- 1
j <- 1

for(i in 1:12){
  ntr <- list.files(pattern = "*wIntrp.tif")
  ntr
	sw <- raster(soilwater[i])
	tm <- raster(temp[i])
	#sw <- aggregate(sw, fact = 5)
	sw.tm <- resample(tm, sw)
	sw.tm <- sw.tm + ((min(val.min.pre)^2)^0.5)
	sw.tm <- sw.tm/temp.range
	sw.tm[sw.tm < 0] <- 0
	sw <- sw/soil.max
	sw.tm <- sw.tm * sw
	for(j in 1:length(ntr)){
		ntr.rast <- raster(ntr[j])
		ntr.rast <- ntr.rast * sw.tm
		bkgrnd <- raster("bkgrnd_2021-06-15.tif")
		ntr.rast <- focal(ntr.rast, w = matrix(1,5,5), mean, na.rm = T)
		ntr.rast <- ntr.rast + bkgrnd
		ntr.name <- gsub('.{4}$','', ntr[i])
		pdf(file = paste(ntr.name, "_Month-", i,".pdf", sep = ""), width = 12, height = 6)
		plot(ntr.rast, col = cols, maxpixels = ncell(ntr.rast))
		dev.off()
		writeRaster(ntr.rast, filename = paste(ntr.name, "_Month-", i,".tif", sep = ""), overwrite = T)
		##### Clear temporary files
		unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
	}
}



