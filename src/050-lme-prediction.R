#------------------------------------------------------------------------------
# Name: 050-lme-prediction
# Author: Chris Reudenbach
# Description:  
# 
#------------------------------------------------------------------------------



# 0 - load packages
#------------------

require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice","rgrass7","link2GI","lubridate","doRNG")

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

forest_structure_metrics <- rast(predstack_forest_metricsFN)  
clim_stat <- readRDS(gw_clim_rdsFN)
clim_stat$doy <- as.numeric(as.character(clim_stat$date_time_hourly, format = "%j"))
clim_stat$hour <- as.numeric(substr(clim_stat$date_time_hourly, 12, 13))

trainDFclean <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))

trainSites <- read_sf(file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))

forest_div <- read_sf(file.path(envrmt$path_auxdata, "forest_divisions.shp"))

# all_lme = readRDS(file.path(envrmt$path_models, "/all_lme.rds"))
model = readRDS(file.path(envrmt$path_models, "/interaction_lme.rds"))



# 3 - run code
#---------------

## --- crop pred stack to MOF MASK ----

# forest_div <- st_transform(forest_div, epsg_number)
# plot(forest_div[1], reset = FALSE)
# plot(trainSites, add = TRUE, pch = 20,col="green")
# pred_stack<- mask(as(forest_structure_metrics,"Raster"), forest_div)
# plot(pred_stack[[1]])
# raster::writeRaster(pred_stack,file.path(envrmt$path_lidar_raster, "mask_pred_forest_structure.tif"))
#
#pred_stack = stack(file.path(envrmt$path_lidar_raster, "mask_pred_forest_structure.tif"))
pred_stack = stack(subset(forest_structure_metrics, c("zpcum2","tree_cluster","aspect")))
names(pred_stack)


## --- 3.1 predict with all models on period ----

# prediction of doy 178 at hour 13
# more or less randomly selected, a day where loggers have very differing measurement
# predict with all developed models and compare result
use_sp()
link2GI::linkGRASS7(as(pred_stack,"Raster"), gisdbase = envrmt$path_data_GRASS, location = "MOF1")



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
res <- foreach(i = seq(5850:8640)+5850, .packages=c("terra","raster","rgrass7","lubridate"),.verbose = TRUE) %dorng% {
  use_sp()
  cat(i, " of ", nrow(mx),"\n")
  all_pred <- pred_mc(i)
}
stopCluster(cl)





