#------------------------------------------------------------------------------
# Name: 070_ex
# Author: Chris Reudenbach
# Description: extract data at 10 defined positions  
# 
#------------------------------------------------------------------------------



# 0 - load packages
#------------------

require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

#root_folder = "~/edu/agis"

# define  additional packages comment if not needed
appendpackagesToLoad = c("stringr","tibble","lme4", "nlme", "CAST", "ggplot2", "lattice","rgrass7","link2GI","lubridate","doRNG","magrittr")

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
source(paste0(root_folder,"/src/000-rspatial-setup.R"))

# extraction radius
ext_rad= 3.0

# create pos data
name = c("Net_1","Net_2","Net_3","Net_4","Net_5","Net_6","Net_7","Net_8","Net_9","Net_10")
Y=c(50.842050,	50.837944,	50.840893,	50.837427,	50.841399,	50.835716,	50.8379259,	50.842949,	50.839332	,50.843793)
X=c(8.685063,	8.678754,	8.674540,	8.673987,	8.670559,	8.680176,	8.6809621,	8.678719,	8.676738,	8.680434)
posDF= data.frame(name,X,Y)
pos_sf = st_as_sf(posDF,coords = c("X", "Y"), crs = "EPSG:4326")
installpos_sf = st_transform(pos_sf,crs = epsg_number)
mapview(pos_sf)

predstack_forest_metricsFN = file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif")
gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")


# - rename files

#mod_type= "_interact1_temp.tif"
mod_type= "_pred_rf_hum.tif"
mod_type= "_pred_rf_temp.tif"
o_year = "2021-01-01" #c("2018-01-01","2019-01-01","2020-01-01","2021-01-01")
year=substr(o_year,1,4)
for (ii in 1){
  temp_files <- list.files(paste0(envrmt$path_data_lev2,year[ii]),pattern = paste0("*",year[ii],".tif"),full.names = T)
  if (length(temp_files) > 0){
    for (i in 1:length(temp_files)){
      split=stringr::str_split(basename(temp_files[i]),pattern = "_",simplify = T)
      hr <-  as.integer(stringr::str_split(split[1,4],pattern = "\\.")[[1]][1])
      dy  = as.integer(split[1,3])
      nd = ymd(as.Date(dy, origin=as.Date(o_year[ii])))
      hour(nd)=hr
      if (hr==0) nd =paste0(nd,"T00:00:00")
      # print(hr)      
      # print(temp_files[i])
      print(paste0(envrmt$path_data_lev2,year[ii],"/",paste0(gsub(nd,pattern = " ",replacement = "T"),mod_type)))
      file.rename(from = temp_files[i]  ,
                  to = paste0(envrmt$path_data_lev2,year[ii],"/",paste0(gsub(nd,pattern = " ",replacement = "T"),mod_type))
      )
    }
  }
  # extract pos timeseries
  extractFN <- list.files(paste0(envrmt$path_data_lev2,year[[ii]]),pattern = "*.tif",full.names = T,recursive = T,include.dirs = T)
  # create monthly files
  for (mo in c("-01-","-02-","-03-","-04-","-05-","-06-","-07-","-08-","-09-","-10-","-11-","-12-")){
    flist = grep(extractFN,pattern = mo)
    month_stack = c(rast(extractFN[flist]))
    ex_values = exactextractr::exact_extract(month_stack, st_buffer(pos_sf, dist = ext_rad),
                                             include_cols="name")
    # DF = dplyr::bind_rows(ex_values) 
    
    extr = dplyr::bind_rows(ex_values) %>% 
      group_by(name) %>% 
      summarise_all(median,na.rm=TRUE)
    
      #rm(month_stack)
    names(extr) = substr(gsub(names(extr), pattern = "temp_", replacement = ""), start = 0, stop = 19)
    #st_write(st_as_sf(ex_values), paste0(envrmt$path_data_lev2,"temp_",mo,".GPKG"), append=FALSE)
    #plot(st_as_sf(ex_values))
    transposed_values = data.frame(t(extr))
    #rm(ex_values)
    transposed_values = data.frame(transposed_values) %>% rownames_to_column('date')
    colnames(transposed_values) =c("datetime",paste0("Net_",seq(10)))
    transposed_values=transposed_values[-1,]
    datetime= ymd_hms(as.character(transposed_values$datetime))
    transposed_values %<>% mutate_if(is.character,as.numeric)
    transposed_values = round(transposed_values[,2:ncol(transposed_values)],digits = 2)
    transposed_values = cbind(datetime = datetime[1:nrow(transposed_values)], transposed_values)
    summary(transposed_values)
    data.table::setcolorder(transposed_values, c("datetime",paste0("Net_",seq(10))))
    write.table( transposed_values, dec = ".", sep = ";",append = FALSE,col.names =  TRUE,,row.names = F,
                 file = paste0(envrmt$path_data_lev2,year[[ii]],"/",substr(mod_type,2,nchar(mod_type)-4),"_",year[[ii]],substr(mo,1,3),".csv"))
    rm(transposed_values)
  }
}

