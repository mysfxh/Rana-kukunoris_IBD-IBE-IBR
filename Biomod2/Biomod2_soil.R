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


############################################################
# 2. Read species occurrence data
############################################################

# Read occurrence records.
# The CSV file should contain three columns:
# Species, X, and Y.
# X = longitude, Y = latitude.
setwd("/home/xiongh/Biomod2/new/result_1//")
points <- read.csv(file = "/home/xiongh/Biomod2/new/result_1/enmeval/GYLW_all.csv", header = T);

# Add a response column.
# Here, all records are presence records, so the response value is 1.
points <- cbind(points, rep.int(1, length(nrow(points)))); 

# Rename columns.
colnames(points) <- c("Species", "X", "Y", "Response");

############################################################
# 3. Read environmental raster layers
############################################################
single_envt <- raster("/home/xiongh/Biomod2/new/result_1/env/saline-alkali soil.asc")

# Read environmental raster layers as a RasterStack.
envt.st <- stack(single_envt); 
myRespName <- 'Xionh'

print("Environmental layers were loaded successfully.")

############################################################
# 4. Format data for biomod2
############################################################

# Format occurrence and environmental data for biomod2.
# Presence records are coded as 1.
# Pseudo-absence points are generated randomly
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

# Set the working directory for model output
setwd("/home/xiongh/Biomod2/new/result_1/")
user_maxent <- list(
  "_allData_allRun" = list(
    path_to_maxent.jar = "./Maxent1",  # Path to the MaxEnt jar file or folder
    memory_allocated = 1024,           # Memory allocated to MaxEnt, in MB
    initial_heap_size = NULL,          # Initial Java heap size
    max_heap_size = NULL,              # Maximum Java heap size
    background_data_dir = "default",   # Background data setting
    visible = FALSE,                   # Do not show the MaxEnt interface
    linear = TRUE,                     # Use linear features
    quadratic = TRUE,                  # Use quadratic features
    product = FALSE,                   # Do not use product features
    threshold = FALSE,                 # Do not use threshold features
    hinge = TRUE,                      # Use hinge features
    betamultiplier = 5,                # Regularization multiplier
    beta_lqp = -1.0,                   # Default regularization for LQP features
    beta_threshold = -1.0,             # Default regularization for threshold features
    beta_hinge = -1.0,                 # Default regularization for hinge features
    defaultprevalence = 0.5            # Default prevalence
  )
)



############################################################
# 6. Create biomod2 modelling options
############################################################

# Pass the user-defined MaxEnt parameters to biomod2.
opt_maxent <- bm_ModelingOptions(
  data.type = 'binary',                      # 数据类型：二进制
  models = c('MAXENT'),                      # 选择 MAXENT 模型
  strategy = 'user.defined',                 # 使用自定义参数
  user.val = list('MAXENT.binary.MAXENT.MAXENT' = user_maxent),  # 使用正确的命名
)

############################################################
# 7. Run MaxEnt modelling
############################################################

# Run MaxEnt models with k-fold cross-validation.
opt_maxent@options$'MAXENT.binary.MAXENT.MAXENT'
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = 'r1_saline-alkali soil',
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
write.csv(eval_df,"./r1saline-alkali soil.csv")

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
  proj.name = 'Maxent_Projection_r1saline-alkali soil',    
  selected.models = "MAXENT",        
  output.format = '.tif',  
  metric.binary = 'TSS',
  metric.filter = 'TSS',
  build.clamping.mask = TRUE,
  models.chosen = max_tss_model
)
