#------------------------------------------------------------------------------
# Name: XXX_RF_model_building
# Author: Lena Perzlmaier
# Description:  
# 
#------------------------------------------------------------------------------

# --- 0 - load packages ----

# ---- project setup ----

require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice","sperrorest","doParallel","CAST")


appendProjectDirList =  c("data/lidar/",
                          "data/GRASS/",
                          "data/data_lev1/models/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/pmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/tmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/topo/",
                          "data/measurements/measure_org",
                          "data/measurements/measure_lev0",
                          "data/measurements/measure_lev1")

# MANDANTORY: calling the setup script also DO NOT change this line
source(file.path(envimaR::alternativeEnvi(root_folder = root_folder),"src/000-rspatial-setup.R"))

trainDF <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))


#trainDF = readRDS(file.path(envrmt$path_auxdata, "trainDFmcScaled.rds"))
# shapefile with tree coordinates
trainSites <- read_sf(file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))

rf_names = c("temp","chm","slope","aspect","tree_cluster","zmax","zskew","zentropy","pzabovezmean","pzabove2","itot","imax","isd","p1th","p2th","p3th","p4th",
"p5th","pground","coverage_fraction","rad_sw_in","rad_sw_out","rad_lw_out","Huma_10m","Huma_2m","Windspeed_10m",
"Wind_direction_10m","rad","cloudiness","Ta_10m","Ta_2m","LAD_metrics_MOF","zpcum2","dtm")


set.seed(100)
predictors = names(trainDF[ , names(trainDF) %in% rf_names])
response <- "rel_humidity"

trainDF = trainDF  %>%
filter(tree_id %in% trainDF$tree_id) %>% 
  mutate(doy_hour_id = paste(doy, hour, tree_id, sep = "_")) 

# --- 3 - start code ----

# partition data spatially with partitions_kmeans
trainSites <- filter(trainSites, trainSites$tree_id %in% trainDF$tree_id)

x_partition = as.vector(sperrorest::partition_kmeans(st_coordinates(trainSites), 
                                         coords = c("X", "Y"), 
                                         nfold = 5, seed1 = 1, 
                                         return_factor = TRUE,  
                                         repetition = 1, 
                                         balancing_steps = 1))[[1]]

partitionDF = data.frame(tree_id = trainSites$tree_id,
                         partition = x_partition)

trainDF = merge(trainDF, partitionDF, by = "tree_id")
trainSites = inner_join(trainSites, partitionDF, by = c("tree_id" = "tree_id"))
mapview(trainSites,zcol="partition")

set.seed(100)
Index <- createDataPartition(trainDF$partition, p = 0.005, list = FALSE)
trainDF_small <- trainDF[Index,]
# create space (10 cluster) time (10 days) folds
trainDF_small$timevar = substr(trainDF_small$doy_hour_id,1,2)
folds <- CreateSpacetimeFolds(trainDF_small, spacevar = "partition", k = 5,timevar = "timevar")
summary(trainDF_small)
cloud <- vector(length = nrow(trainDF_small))
for (i in 1:nrow(trainDF_small)) {
  if(trainDF_small[i,]$rad <= 0){
    cloud[i] <- 0
  } else {
    cloud[i] <- trainDF_small[i,]$rad_sw_in/trainDF_small[i,]$rad
    if (cloud[i] <=0) cloud[i] <- 0
  }
}

##--- filter zero or near-zero values
#  nzv = nearZeroVar(subset(trainDF_small, select=-c(cst_id,coverage_fraction,tree_id,doy_hour_id,day,partition,timevar,doy,hour)))
#  names(trainDF_small[,nzv])
# # if (length(nzv) > 0) traintmp = traintmp[, -nzv]
# 
# ##--- filter correlations that are > cor_cutoff
# filt = findCorrelation(cor(subset(trainDF_small, select=-c(cst_id,coverage_fraction,tree_id,doy_hour_id,day,partition,timevar,doy,hour)), use = "complete"), cutoff = 0.98)
# #traintmp = traintmp[,-filt]


trainDF_small$cloudiness <- cloud
# forward feature selection
cl <- makeCluster(8)
registerDoParallel(cl)

ctrl <- trainControl(method = "cv", 
                     index = folds$index,
                     indexOut = folds$indexOut,
                     savePredictions = TRUE)

ffs_model <- CAST::ffs(trainDF_small[,predictors],
                       trainDF_small[,response],
                       method = "rf",
                       metric="RMSE",
                       importance =TRUE,
                       tuneGrid = expand.grid(mtry = 2),
                       #ntree = 100,
                       trControl = ctrl)
stopCluster(cl)

ffs_model$selectedvars


saveRDS(ffs_model, file.path(envrmt$path_auxdata, "rf_ffs_base_hum.rds"))
ffs_model = readRDS(file.path(envrmt$path_auxdata, "rf_ffs_base_hum.rds"))

# 2 - final hyperparameter tuning

predictors = ffs_model$selectedvars

ctrl <- trainControl(method = "cv", 
                     savePredictions = TRUE,
                     index = folds$index, 
                     indexOut = folds$indexOut)

tunegrid = expand.grid(.mtry = c(2:4))

modellist <- list()
cl <- makeCluster(30)
registerDoParallel(cl)

for (ntree in c(200, 400, 600)) {
  
  set.seed(100)

  fit = train(trainDF_small[,predictors],
              trainDF_small[,response],
              method = "rf",
              metric = "RMSE",
              ntree = ntree,
              tuneGrid = tunegrid,
              trControl = ctrl)
 
   key <- toString(ntree)
  
  modellist[[key]] <- fit
}

stopCluster(cl)

rf_tuned = resamples(modellist)
summary(rf_tuned)
plot(rf_tuned)

# best mtry is 5
# highest accuracy was calculated with ntree = 600

saveRDS(rf_tuned, file.path(envrmt$path_auxdata, "rf_ffs_base_hum_final_tuning.rds"))
rf_tuned = readRDS(file.path(envrmt$path_auxdata,"rf_ffs_base_hum_final_tuning.rds"))

# 4 - final model

predictors = rf_ffs$selectedvars

ctrl <- trainControl(method = "cv", 
                     savePredictions = TRUE,
                     index = folds$index, 
                     indexOut = folds$indexOut)

mtry = 3
tunegrid = expand.grid(.mtry = mtry)

ntree = 600

cl = makeCluster(16)

registerDoParallel(cl)
  
set.seed(100)
  
rf_final = train(trainDF_small[,predictors],
              trainDF_small[,response],
              method = "rf",
              metric = "RMSE",
              ntree = ntree,
              tuneGrid = tunegrid,
              trControl = ctrl)


stopCluster(cl)

rf_final

plot(varImp(rf_final))

saveRDS(rf_final, file.path(envrmt$path_auxdata, "rf_ffs_base_hum_final.rds"))
rf_final = readRDS(file.path(envrmt$path_auxdata,"rf_ffs_base_hum_final.rds"))
