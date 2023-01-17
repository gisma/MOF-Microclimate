#------------------------------------------------------------------------------
# Name: 040-make-RF-model
# Author: Chris Reudenbach
# Description:  training and tuning of rf model with ffs and ltllo out 
# 
#------------------------------------------------------------------------------

# ---- project setup ----

require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()
root_folder="~/edu/agis"
# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice","sperrorest","doParallel","CAST","lubridate")


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

# controlling the run settings either "temp" or "rel_humidity"
response <- "temp"
rf_names = c("chm","slope","aspect","tree_cluster","zmax","zentropy","pzabovezmean","pzabove2"
             ,"rad","cloudiness","Ta_10m","VCI_metrics_MOF","entropy_metrics_MOF","LAD_metrics_MOF" ,"zpcum2","dsm")

trainDF <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))


#trainDF = readRDS(file.path(envrmt$path_auxdata, "trainDFmcScaled.rds"))
# shapefile with tree coordinates

trainSites <- sf::read_sf(file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))

set.seed(100)
predictors = names(trainDF[ , names(trainDF) %in% rf_names])


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
                                                     repetition =1, 
                                                     balancing_steps = 1))[[1]]

partitionDF = data.frame(tree_id = trainSites$tree_id,
                         partition = x_partition)

trainDF = merge(trainDF, partitionDF, by = "tree_id")
trainSites = inner_join(trainSites, partitionDF, by = c("tree_id" = "tree_id"))
mapview(trainSites,zcol="partition")


set.seed(100)
Index <- createDataPartition(trainDF$partition, p = 0.005, list = FALSE)
trainDF_small <- trainDF[Index,]
trainDF_small_test <- trainDF[-Index,]
# create space (10 cluster) time (10 days) folds
trainDF_small$timevar = substr(trainDF_small$doy_hour_id,1,3)
folds <- CreateSpacetimeFolds(trainDF_small, spacevar = "partition", timevar = "timevar",k=5)
summary(trainDF_small[,predictors])
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
# 
trainDF_small$cloudiness <- cloud


# forward feature selection first run with ntree = 100
cl <- makePSOCKcluster(16)
registerDoParallel(cl)

ctrl <- trainControl(method = "cv", 
                     index = folds$index,
                     indexOut = folds$indexOut,
                     savePredictions = TRUE)
tunegrid = expand.grid(.mtry = c(2:floor(sqrt(length(predictors)))))
ffs_model <- CAST::ffs(trainDF_small[,predictors],
                       trainDF_small[,response],
                       method = "rf",
                       metric="RMSE",
                       importance =TRUE,
                       tuneGrid = tunegrid,
                      # ntree = 75,
                       trControl = ctrl,
                       allowParallel=TRUE)
stopCluster(cl)

ffs_model$selectedvars
plot(varImp(ffs_model))

saveRDS(ffs_model, paste0(envrmt$path_auxdata, "rf_ffs_base_",response ,".rds"))
ffs_model = readRDS(paste0(envrmt$path_auxdata, "rf_ffs_base_",response ,".rds"))



#-- second run parameter tuning
predictors = ffs_model$selectedvars
ctrl <- trainControl(method = "cv", 
                     savePredictions = TRUE,
                     index = folds$index, 
                     indexOut = folds$indexOut)
tunegrid = expand.grid(.mtry = c(2:floor(sqrt(length(predictors)))))
modellist <- list()

cl <- makeCluster(30)
registerDoParallel(cl)

for (ntree in seq(from=200,to=800,by=100)) {
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
plot(varImp(rf_tuned))

saveRDS(rf_tuned, paste0(envrmt$path_auxdata, "rf_ffs_base_",response ,"_final_tuning.rds"))
rf_tuned = paste0(file.path(envrmt$path_auxdata,"rf_ffs_base_",response ,"_final_tuning.rds"))

## get best mtry and ntree
# for humidity: mtry = 2 ; ntree = 475
# for temp mtry = 2, ntree =550
if (response == "temp"){
  mtry = 2
  ntree = 550
  } else if (response == "rel_humidity"){
  mtry = 2
  ntree = 475
}

# train final model
predictors = rf_ffs$selectedvars
ctrl <- trainControl(method = "cv", 
                     savePredictions = TRUE,
                     index = folds$index, 
                     indexOut = folds$indexOut)
tunegrid = expand.grid(.mtry = mtry)

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

plot(varImp(rf_final))
rf_final$l
saveRDS(rf_final, paste0(envrmt$path_auxdata, "rf_ffs_base_",response ,"_final.rds"))
rf_final = readRDS(paste0(envrmt$path_auxdata,"rf_ffs_base_",response ,"_final.rds"))

#rmse <- round(sqrt(mean((cv$pred - cv$obs)^2, na.rm = TRUE)), 2)

# validation of the model
# get all cross-validated predictions and calculate kappa
