############################################################
# MaxEnt species distribution modelling using biomod2
############################################################


library(terra)
library(biomod2)
library(raster)
library(gam)
library(tidyterra)
library(ggplot2)
library(dismo)
library(doParallel)
#setwd("D://BaiduSyncdisk//enmevalu7te_for_enmeval_2.0//enmevaluate_for_enmeval_2.0//")
#setwd("/home/xiongh/Biomod2/")
setwd("/home/xiongh/Biomod2/new/result_1//")
points <- read.csv(file = "/home/xiongh/Biomod2/new/result_1/enmeval/GYLW_all.csv", header = T);
points <- cbind(points, rep.int(1, length(nrow(points))));
colnames(points) <- c("Species", "X", "Y", "Response");

single_envt <- raster("/home/xiongh/Biomod2/new/result_1/env/aspect.asc")

envt.st <- stack(single_envt); 
myRespName <- 'Xionh'


print("yes")

myBiomodData <- BIOMOD_FormatingData(
  resp.var = ifelse(points[, 4] == 1, 1, NA),  
  resp.xy = points[, 2:3],
  resp.name = as.character(points[1, 1]),
  expl.var = envt.st,
  PA.nb.rep = 10,         
  PA.nb.absences = 1000,    
  PA.strategy = "random",   
  filter.raster = TRUE

)



# 定义 MAXENT 模型的自定义参数
#setwd('D://Maxent//')
setwd("/home/xiongh/Biomod2/new/result_1/")
user_maxent <- list(
  '_allData_allRun' = list(
    path_to_maxent.jar = './Maxent1', 
    memory_allocated =1024,                    
    initial_heap_size = NULL,             
    max_heap_size = NULL,                 
    background_data_dir = "default",       
    visible = FALSE,                            
    linear = TRUE,                              
    quadratic = TRUE,                           
    product = FALSE,                             
    threshold = FALSE,                           
    hinge = TRUE,                              
    betamultiplier = 5,                         
    beta_lqp = -1.0,                            
    beta_threshold = -1.0,                      
    beta_hinge = -1.0,                          
    defaultprevalence = 0.5                     
  )
)



opt_maxent <- bm_ModelingOptions(
  data.type = 'binary',                      
  models = c('MAXENT'),                      
  strategy = 'user.defined',                 
  user.val = list('MAXENT.binary.MAXENT.MAXENT' = user_maxent),  
)


opt_maxent@options$'MAXENT.binary.MAXENT.MAXENT'
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = 're1_75apect',
  models = c('MAXENT'),
  CV.strategy = 'kfold',
  CV.nb.rep = 10,
  CV.perc = 0.75,
  CV.k = 5,
  OPT.user =opt_maxent,  
  metric.eval = c('TSS', 'ROC'),
  do.full.models = FALSE,
  nb.cpu = 200
)

eval_df <- get_evaluations(myBiomodModelOut)
write.csv(eval_df,"./re1_aspect.csv")
#eval_df<-get_evaluations(myBiomodModelOut)
tss_values <- eval_df[eval_df$metric.eval == "TSS", ]
print(tss_values[, c("run", "calibration", "validation", "evaluation")])
max_tss_row <- tss_values[which.max(tss_values$validation), ]
max_tss_model <- max_tss_row$full.name
print(max_tss_model)

my_projection <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, 
  new.env = envt.st ,                  
  proj.name = 'Maxent_Projection_re1aspect',     
  selected.models = "MAXENT",          
  output.format = '.tif',  
  metric.binary = 'TSS',
  metric.filter = 'TSS',
  build.clamping.mask = TRUE,
  models.chosen = max_tss_model
)
