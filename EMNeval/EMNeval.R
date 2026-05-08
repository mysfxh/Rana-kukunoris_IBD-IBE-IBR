library(ENMeval)
library(raster)
library(rJava)
library(dismo)
library(magrittr)
library(sp)
library(doParallel)
options(java.parameters = "-Xmx64g")
options(expressions = 500000)

env.files <- list.files(path ="/home/xiongh/Biomod2/new/result_1/env/", pattern = ".asc", full.names = TRUE)
env <- stack(env.files)

occ <- read.csv("/home/xiongh/Biomod2/new/result_1/enmeval/GYLW_all.csv")[,-1]
colnames(occ)<-c("x","y")
setwd("/home/xiongh/Biomod2//Maxent/")
length(which(!is.na(values(subset(env, 1)))))
sapply(1:nlayers(env), function(i) sum(is.na(values(env[[i]]))))
bg <-dismo::randomPoints(env[[1]],n=10000)%>%as.data.frame()
enmeval_results2 <- ENMevaluate(occ, env, bg = bg[,1:2], tune.args = list(fc = c("L","Q","H","LQ","LQH", "LQHP","LQHPT"), rm =  c(0.1,seq(0.5,5,0.5))), partitions = "randomkfold",algorithm = "maxent.jar",categoricals = names(env)[c(10)],clamp = TRUE)

write.csv(enmeval_results2@results,"/home/xiongh/Biomod2/new/result_1/enmeval/enm_2914new.csv")
