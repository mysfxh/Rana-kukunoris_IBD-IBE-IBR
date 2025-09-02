library(ENMTools)
library(raster)
library(raster)
library(sp)
library(terra)

setwd("/home/xiongh/Biomod2/env/")
envtList <- list.files(pattern = ".asc");
envt.st <- stack(envtList);
cor <- raster.cor.matrix(envt.st,method = "pearson")
write.csv(cor,"./cor.csv")

