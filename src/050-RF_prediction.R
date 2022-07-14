require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice","sperrorest","doParallel","CAST","doRNG")

# define additional subfolders comment if not needed
  
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

predstack_forest_metricsFN = file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif")
gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")



# 2 load files
#--------------

model_rf <- readRDS(file.path(envrmt$path_auxdata,"rf_ffs_base_hum_final.rds"))

forest_structure_metrics <- rast(predstack_forest_metricsFN)  
clim_stat <- readRDS(gw_clim_rdsFN)
clim_stat$doy <- as.numeric(as.character(clim_stat$date_time_hourly, format = "%j"))
clim_stat$hour <- as.numeric(substr(clim_stat$date_time_hourly, 12, 13))

trainDFclean <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))

trainSites <- read_sf(file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))

forest_div <- read_sf(file.path(envrmt$path_auxdata, "forest_divisions.shp"))

pred_stack = stack(subset(forest_structure_metrics, c("zskew","p4th","zentropy","zpcum2")))
names(pred_stack)
use_sp()
link2GI::linkGRASS7(as(pred_stack,"Raster"), gisdbase = envrmt$path_data_GRASS, location = "MOF1")


# 3 - run code
#--------------


# perform some predictions
# make matrix containing 365*24 doy and hours with corresponding dates for the prediction period
ts = c(as.POSIXct("2099-01-01 00:00:00 CET"),seq.POSIXt(min(clim_stat$date_time_hourly,na.rm = TRUE),max(clim_stat$date_time_hourly,na.rm = TRUE),by="hour"))
i=1
j=1
k=1
mm = list()
for (h in c(2018,2019,2020,2021)){
  ndays = 365
  if (h == 2020) ndays = 366
  for (j in (seq(ndays))){
    for (k in seq(1:24)){
      if (j>365) j= 365  
      mm[[i]] =  data.frame(date= ts[i+1], day=j, hour=k-1)
      i=i+1
      k=k+1
    }
  }
}
mx =do.call(rbind,mm)


# predict
cl <- makeCluster(28)
registerDoParallel(cl)
kind <- RNGkind()
set.seed(1234, kind = kind[1])
res <- foreach(i = seq(1:8640), .packages=c("terra","raster","rgrass7","lubridate","caret","randomForest"),.verbose = TRUE) %dorng% {
  use_sp()
  cat(i, " of ", nrow(mx),"\n")
  all_pred <- pred_rf(i)
}
stopCluster(cl)
