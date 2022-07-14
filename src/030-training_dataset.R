#
# Name: 30_training_dataset
# Author: Lena Perzlmaier, Chris Reudenbach
# Description:  creates the training data frame
#             - stack topo and forest structure  -> pred_stack
#             - get coordinates of the trainSites and extract information of predStack
#             - merge information of the step before with data from the climate station
#             - merge with radiation information
# Input:   see source files
             # Output: - pred_stack "C:/Users/Lena/Documents/edu/mpg-envinsys-plygrnd/data/predictors/mc_pred_stack.tif"
#         - trainDF "C:/Users/Lena/Documents/edu/mpg-envinsys-plygrnd/data/auxdata/trainDFmc.rds"
# 
#

# ---- project setup ----
# ----project setup----
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lidR","future","lwgeom","tmap","ggplot2","VIM","Hmisc","mice")

# define additional subfolders comment if not needed
tmpPath = paste0("data/lidar/",prefix,"/") 
appendProjectDirList =  c("data/lidar/",
                          "data/GRASS/",
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

gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")
gw_clim_csvFN = file.path(envrmt$path_measure_org, "302_KlimaturmCaldernWiese_complete.csv")
combo_clim_rdsFN = file.path(envrmt$path_measure_lev0, "climate_stations_combined.rds")
# ---- read preprocessed data ----
# shapefile with tree coordinates
trainSites <- read_sf(file.path(envrmt$path_auxdata, "core_study_trees.shp"))
trainSites = st_transform(trainSites,crs = epsg_number)

# predstack
#topo <- stack(file.path(envrmt$path_dem,"pred_topo.tif"))
#forest <- readRDS(file.path(envrmt$path_MOF_lidar_2018,"pred_forest_structure.rds"))
pred_stack = rast(file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif"))
# climate station Grubenwiese
# trainDF
trainDF <- readRDS(combo_clim_rdsFN)
head(trainDF)
tail(trainDF)

clim_stat <- readRDS(gw_clim_rdsFN)

# location coordinate of the climate station

cstation = st_sfc(st_point(cbind(8.6832, 50.8405)),crs = 4326)
cstation = st_transform(st_sf(data.frame(tree_id="grubenwiese", geom=cstation)),crs = 25832)
 t <-trainSites %>% 
   st_as_sf(coords = c("Easting", "Northing"), crs = 25832L, remove = FALSE)
 t= trainSites[1,]
 t$id = 999; t$tree_id = "gw_station";t$species_id="NA";t$northing=5632136;t$easting=477694.2;t$geometry=cstation$geometry;t$color_rgb="0,0,255";t$remarks="campbell";t$state_rema="campbell"
 trainSites=rbind(trainSites,t)
 st_write(trainSites,file.path(envrmt$path_auxdata, "core_study_trees_plus.shp"))
# ---- 2  run code ----

# ---- prepare trainDF ----
# day 
trainDF$doy <- as.numeric(as.character(trainDF$date, format = "%j"))
# hour
trainDF$hour <- as.numeric(substr(trainDF$date, 12, 13))

# ---- 2.3 - prepare trainSites ----
# relevant tree IDs  GW climate station will be merged later to EACH TT
IDs <- unique(c(trainDF$tree_id,"gw_station"))

# filter shapefile to only containt relevant tree IDs
trainSites <- filter(trainSites, trainSites$tree_id %in% IDs)
trainSites <- st_crop(trainSites, st_bbox(pred_stack))
mapview(as(pred_stack[[1]],"Raster"),fgb=F)+ trainSites

# ---- 2.4 - extract data ----
extr = exactextractr::exact_extract(pred_stack, st_buffer(trainSites,dist = 3),  force_df = TRUE,
                                    include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "tree_id") 

extr = dplyr::bind_rows(extr) %>%
  group_by(tree_id) %>% 
  summarise_all(median,na.rm=TRUE)


trainDF2 <- full_join(trainDF, extr)

# less observation than before because cst_id with NAs where removed

# ---- 2.5 - merge with climate station data ----
#trainDF <- right_join(trainDF2, clim_stat, by.x = "date", by.y = "date_time_hourly")
trainDF <- full_join(trainDF2, clim_stat, by = c("date" = "date_time_hourly"))
# add label col
trainDF <- trainDF %>% mutate(doy_hour_id = paste(doy, hour, tree_id,  sep = "_"))
#colnames(trainDF)[which(names(trainDF) == "cst_id.x")] = "cst_id"

ggplot(data = trainDF, aes(x = date, y = Ta_2m)) +
  geom_point(size = 0.25)+
  facet_wrap(~tree_id)

# ---- 2.6 get GRASS radiation information ----

use_sp()
link2GI::linkGRASS7(as(pred_stack,"Raster"), gisdbase = envrmt$path_data_GRASS, location = "MOF1")

doy_hour <- unique(trainDF[,c("doy", "hour")])
doy_hour$hour_chr <- sprintf("%02d", doy_hour$hour)


md.pattern(trainDF)
mice_plot <- aggr(trainDF, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(trainDF), cex.axis=.3
                  ,gap=3, ylab=c("Missing data","Pattern"))

# extract RAD for treetalker positions
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

rad <- foreach(i = seq(nrow(doy_hour)-1),.packages=c("raster","rgrass7","dplyr","sf","exactextractr"),.verbose = TRUE) %dopar% {
  use_sp()
  file_name <- paste("tot_rad_", doy_hour[i,1], "_", doy_hour[i,3], ".00", sep = "")
  ras <- raster(rgrass7::read_RAST(file_name))
  cst_ids <- filter(trainDF, doy == doy_hour[i,1] & hour == doy_hour[i,2])
  cst_ids <- c(unique(cst_ids$tree_id),"gw_station")
  trainSites_filtered <- filter(trainSites, tree_id %in% cst_ids)
  if (nrow(trainSites_filtered) > 0){
  extr = exactextractr::exact_extract(ras, st_transform(st_buffer(trainSites_filtered,dist = 3.),crs = st_crs(ras)),
                                      force_df = TRUE,
                                      include_cols = "tree_id",
                                      include_area =T) 
  extr = dplyr::bind_rows(extr) %>%
    group_by(tree_id) %>% 
    summarise(mean = mean(value), n = n())
  
  extr <- data.frame(tree_id = extr$tree_id,
                     doy = doy_hour[i,1],
                     hour = doy_hour[i,2],
                     rad = extr$mean)
  }
}

stopCluster(cl)
rad <- do.call("rbind", rad)


# merge results with trainDF
trainDF2 <- merge(trainDF, rad, by = c("tree_id", "doy", "hour"))
plot(trainDF2$rad, trainDF2$temp)
plot(trainDF2$rad, trainDF$Ta_2m)

# ---- save training data frame ----
saveRDS(trainDF2, file.path(envrmt$path_measure_lev0, "trainDFmc.rds"))
trainDF <- readRDS(file.path(envrmt$path_measure_lev0, "trainDFmc.rds"))
