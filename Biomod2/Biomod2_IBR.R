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
points <- read.csv(file = "/home/xiongh/Biomod2/GYLW_all.csv", header = T);
points <- cbind(points, rep.int(1, length(nrow(points)))); #新生成一列代表当前存在的点
colnames(points) <- c("Species", "X", "Y", "Response");

# 假设你有环境变量数据 (如 BIOCLIM 数据)，存储在 myExpl 中
envtList <- list.files(path = "/home/xiongh/Biomod2/new/result_1/env/IBR/",pattern = ".asc");
print(envtList)
setwd("/home/xiongh/Biomod2/new/result_1/env/")
envt.st <- stack(envtList); # 读取栅格数据
# 格式化数据
print("yes")

myBiomodData <- BIOMOD_FormatingData(
  resp.var = ifelse(points[, 4] == 1, 1, NA),  # 将出现点的值设置为1，其他点设置为NA
  resp.xy = points[, 2:3],
  resp.name = as.character(points[1, 1]),
  expl.var = envt.st,
  PA.nb.rep = 10,           # 设置生成伪缺失数据的重复次数
  PA.nb.absences = 1000,    # 设置每次选择1000个伪缺失数据点
  PA.strategy = "random",    # 设置伪缺失数据的选择策略，使用随机选择
  filter.raster = TRUE

)



# 定义 MAXENT 模型的自定义参数
#setwd('D://Maxent//')
setwd("/home/xiongh/Biomod2/new/result_1/")
user_maxent <- list(
  '_allData_allRun' = list(
    path_to_maxent.jar = './Maxent1',  # MaxEnt JAR 文件路径
    memory_allocated =1024,                     # 分配的内存（MB）
    initial_heap_size = NULL,             # 初始堆内存大小
    max_heap_size = NULL,                 # 最大堆内存大小
    background_data_dir = "default",       # 背景数据目录
    visible = FALSE,                            # 是否显示 MaxEnt 界面
    linear = TRUE,                              # 是否使用线性特征
    quadratic = TRUE,                           # 是否使用二次特征
    product = FALSE,                             # 是否使用乘积特征
    threshold = FALSE,                           # 是否使用阈值特征
    hinge = TRUE,                               # 是否使用铰链特征
    betamultiplier = 5,                         # 规则化系数倍数
    beta_lqp = -1.0,                            # 线性、二次和乘积特征的正则化系数
    beta_threshold = -1.0,                      # 阈值特征的正则化系数
    beta_hinge = -1.0,                          # 铰链特征的正则化系数
    defaultprevalence = 0.5                     # 物种的默认发生概率
  )
)


# 将这些自定义参数传入 bm_ModelingOptions
opt_maxent <- bm_ModelingOptions(
  data.type = 'binary',                      # 数据类型：二进制
  models = c('MAXENT'),                      # 选择 MAXENT 模型
  strategy = 'user.defined',                 # 使用自定义参数
  user.val = list('MAXENT.binary.MAXENT.MAXENT' = user_maxent),  # 使用正确的命名
)

# 查看 MAXENT 模型的自定义参数
opt_maxent@options$'MAXENT.binary.MAXENT.MAXENT'
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = 'r1_ibr',
  models = c('MAXENT'),
  CV.strategy = 'kfold',
  CV.nb.rep = 10,
  CV.perc = 0.75,
  CV.k = 5,
  OPT.user =opt_maxent,  # 选择我们自己定义的模型
  metric.eval = c('TSS', 'ROC'),
  do.full.models = FALSE,
  nb.cpu = 200
)
#目前出现一个问题：
eval_df <- get_evaluations(myBiomodModelOut)
write.csv(eval_df,"./r1__ibr.csv")
#eval_df<-get_evaluations(myBiomodModelOut)
tss_values <- eval_df[eval_df$metric.eval == "TSS", ]
print(tss_values[, c("run", "calibration", "validation", "evaluation")])
max_tss_row <- tss_values[which.max(tss_values$validation), ]
max_tss_model <- max_tss_row$full.name
print(max_tss_model)

my_projection <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,  # 模型输出对象
  new.env = envt.st ,                  # 投影的环境数据（RasterStack或RasterBrick格式）
  proj.name = 'Maxent_Projection_r1ibr',     # 项目名称
  selected.models = "MAXENT",          # 选择模型
  output.format = '.tif',  # 指定输出格式为 .asc
  metric.binary = 'TSS',
  metric.filter = 'TSS',
  build.clamping.mask = TRUE,
  models.chosen = max_tss_model
)

