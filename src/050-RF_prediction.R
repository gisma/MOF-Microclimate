require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice","sperrorest","doParallel","CAST","doRNG","lubridate","rgrass7")

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

response = "temp"
temper = TRUE
predstack_forest_metricsFN = file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif")
gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")



# 2 load files
#--------------
model= readRDS(paste0(envrmt$path_auxdata,"rf_ffs_base_temp_final.rds"))
#model_rf <- readRDS(paste0(envrmt$path_auxdata,"rf_ffs_base_",response ,"_final.rds"))
#ffs_model = readRDS(paste0(envrmt$path_auxdata, "rf_ffs_base_hum.rds"))
forest_structure_metrics <- rast(predstack_forest_metricsFN)  
clim_stat <- readRDS(gw_clim_rdsFN)
clim_stat$doy <- as.numeric(as.character(clim_stat$date_time_hourly, format = "%j"))
clim_stat$hour <- as.numeric(substr(clim_stat$date_time_hourly, 12, 13))

trainDFclean <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))

trainSites <- read_sf(file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))

forest_div <- read_sf(file.path(envrmt$path_auxdata, "forest_divisions.shp"))
names(forest_structure_metrics)[1] ="dtm"

# perform some predictions
# make matrix containing 365*24 doy and hours with corresponding dates for the prediction period
ts = c(as.POSIXct("2018-10-08 00:00:00 CET"),seq.POSIXt(as.POSIXct("2021-12-06 22:00:00 CET"),max(clim_stat$date_time_hourly,na.rm = TRUE),by="hour"))
ts= ts[2:length(ts)]

jd = yday(ts)
hr = hour(ts)
mx=data.frame(date=ts,day=jd,hour=hr)

saveRDS(mx,paste0(envrmt$path_auxdata,"mx.rds"))
mx = readRDS(paste0(envrmt$path_auxdata,"mx.rds"))
# 3 - run code
#--------------

#predstack_names= (! predictors %in% c("rad","rad_sw_in","rad_sw_out","rad_lw_out","Huma_10m","Huma_2m","Windspeed_10m","Wind_direction_10m","cloudiness","Ta_10m","Ta_2m"))
predstack_names =c("aspect")

pred_stack = stack(subset(forest_structure_metrics,predstack_names))
names(pred_stack)


if (!temper){
  # predict
  cl <- makeCluster(28)
  registerDoParallel(cl)
  #kind <- RNGkind()
  set.seed(1234)#, kind = kind[1])
  res <- foreach(i = seq(1:nrow(mx)), .packages=c("terra","raster","lubridate","caret","randomForest"),.verbose = TRUE) %dopar% {
    cat(i, " of ", nrow(mx),"\n")
    const <- data.frame(
      Ta_10m = clim_stat$Ta_10m[i],
      Windspeed_10m = clim_stat$Windspeed_10m[i]
    )
    hour_chr <-  sprintf("%02d", mx[i,3])
    file_name <- paste("tot_rad_", mx[i,2], "_", hour_chr, ".00", sep = "")
    yr = lubridate::year(as.POSIXct(mx[i,1]))
    
    pred <- raster::predict(object = pred_stack, 
                            model = model_rf , 
                            allow.new.levels=TRUE, 
                            const=const)
    terra::writeRaster(pred,
                       paste0(envrmt$path_data_lev2,file_name,"_pred_rf_hum_",yr,".tif"),
                       overwrite=TRUE,gdal=c("COMPRESS=DEFLATE"))  
  }
  stopCluster(cl)
} else if (temper) 
{
  # setup GRASS
  use_sp()
  link2GI::linkGRASS7(as(pred_stack,"Raster"), gisdbase = envrmt$path_data_GRASS, location = "MOF1")
  cl <- makeCluster(28)
  registerDoParallel(cl)
  kind <- RNGkind()
  set.seed(1234, kind = kind[1])
  res <- foreach(i = seq(1:nrow(mx)), .packages=c("terra","raster","rgrass7","lubridate","caret","randomForest"),.verbose = TRUE) %dorng% {
    use_sp()
    cat(i, " of ", nrow(mx),"\n")
    all_pred <- pred_rf_temp(i)
  }
  stopCluster(cl)
}
