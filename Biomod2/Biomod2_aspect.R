############################################################
# MaxEnt species distribution modelling using biomod2
############################################################

############################################################
# 1. Load required packages
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

############################################################
# 2. Read species occurrence data
############################################################
points <- read.csv(file = "/home/xiongh/Biomod2/new/result_1/enmeval/GYLW_all.csv", header = T);

# Add a response column.
# Here, all records are presence records, so the response value is 1.
points <- cbind(points, rep.int(1, length(nrow(points))));

# Rename columns.
colnames(points) <- c("Species", "X", "Y", "Response");

############################################################
# 3. Read environmental raster layers
############################################################

single_envt <- raster("/home/xiongh/Biomod2/new/result_1/env/aspect.asc")
# Read environmental raster layers as a RasterStack.
envt.st <- stack(single_envt); 
myRespName <- 'Xionh'


print("Environmental layers were loaded successfully.")


############################################################
# 4. Format data for biomod2
############################################################

# Format occurrence and environmental data for biomod2.
# Presence records are coded as 1.
# Pseudo-absence points are generated randomly.

myBiomodData <- BIOMOD_FormatingData(
  resp.var = ifelse(points[, 4] == 1, 1, NA),  
  resp.xy = points[, 2:3],
  resp.name = as.character(points[1, 1]),
  expl.var = envt.st,
  PA.nb.rep = 10,          # Number of pseudo-absence replicates   
  PA.nb.absences = 1000,   # Number of pseudo-absence points per replicate  
  PA.strategy = "random",  # Random pseudo-absence selection 
  filter.raster = TRUE     # Remove duplicated records in the same raster cell

)



############################################################
# 5. Set custom MaxEnt parameters
############################################################

# Set the working directory for model output.
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

############################################################
# 6. Create biomod2 modelling options
############################################################

# Pass the user-defined MaxEnt parameters to biomod2.

opt_maxent <- bm_ModelingOptions(
  data.type = 'binary',                      
  models = c('MAXENT'),                      
  strategy = 'user.defined',                 
  user.val = list('MAXENT.binary.MAXENT.MAXENT' = user_maxent),  
)

# Check the MaxEnt options.
opt_maxent@options$'MAXENT.binary.MAXENT.MAXENT'

############################################################
# 7. Run MaxEnt modelling
############################################################

# Run MaxEnt models with k-fold cross-validation.
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

############################################################
# 8. Extract and save model evaluation results
############################################################

# Extract model evaluation results.
eval_df <- get_evaluations(myBiomodModelOut)

# Save evaluation results.
write.csv(eval_df,"./re1_aspect.csv")

# Extract TSS values.
tss_values <- eval_df[eval_df$metric.eval == "TSS", ]

# Print TSS results for each model run.
print(tss_values[, c("run", "calibration", "validation", "evaluation")])

# Select the model with the highest validation TSS.
max_tss_row <- tss_values[which.max(tss_values$validation), ]
max_tss_model <- max_tss_row$full.name
print(max_tss_model)

############################################################
# 9. Project the best MaxEnt model
############################################################

# Project the selected MaxEnt model to the environmental layers.
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
