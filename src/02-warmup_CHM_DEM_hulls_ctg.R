#------------------------------------------------------------------------------
# Name: 20_CHM_DEM_hulls_tiles.R
# Type: control script 
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  create CHM, DSM, tree segmentation and metrics
# Data: regular las LiDAR data sets 
# Copyright:GPL (>= 3) 
# Date: 2021-12-10 
# V-2022-04-02; adaption to lidR 4.0 , reorganisation
#------------------------------------------------------------------------------
# 0 - project setup
#-----------------------------
# 0 - specific setup
#-----------------------------
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()
# workaround subfolder
appendpackagesToLoad = c("lidR","future","lwgeom","tmap")
appendProjectDirList =  c("data/lidar/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/")

source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)

#--- further variables
# (re)-calculate: calculate = TRUE
# read data: calculate = TRUE
calculate = TRUE

# filenames 
# NOTE FN are corresponding with the lidrc catalog path settings
ctgFN = paste0(envrmt$path_lidar_ctg,"ctg.rds")
demFN = paste0(envrmt$path_lidar_raster,"rasterize_terrain.vrt")
chmFN = paste0(envrmt$path_lidar_raster,"rasterize_canopy.vrt")
treesFN = paste0(envrmt$path_lidar_vector,"segmentation_sapflow.shp")
sapflow_metricsFN = paste0(envrmt$path_lidar_raster,"segmentation_sapflow.tif")
las_fileFN = paste0(envrmt$path_lidar_orig,"MOF_lidar_2018.las")



# resolution in meter
res = 2.0

#- parallelisation in this case 8 Cores
future::plan(multisession, workers = 12)
set_lidr_threads(12)

#------------------------------------------------------------------------------
# 1 - start processing
#-----------------------------
if (calculate){
  #---- define lidR catalog ----
  #- general catalog settings
  ctg <- lidR::readLAScatalog(las_fileFN)
  projection(ctg) <- 25832
  lidR::opt_chunk_size(ctg) = 750
  lidR::opt_chunk_buffer(ctg) <- 25
  lidR::opt_progress(ctg) <- TRUE
  lidR::opt_laz_compression(ctg) <- TRUE
  ctg@output_options$drivers$Raster$param$overwrite <- TRUE
  ctg@output_options$drivers$SpatRaster$param$overwrite <- TRUE
  ctg@output_options$drivers$SpatVector$param$overwrite <- TRUE
  ctg@output_options$drivers$stars$param$overwrite <- TRUE
  ctg@output_options$drivers$sf$param$overwrite <- TRUE
  #---- height normalisation (point cloud) ----
  # source: https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
  # first derive dem
  lidR::opt_output_files(ctg) <- paste0(envrmt$path_lidar_raster,"dem/{ID}","_dem")
  dem <- rasterize_terrain(ctg, res = res,lidR::knnidw(k = 6L, p = 2))
  
  #- second normalize height (point cloud)
  lidR::opt_output_files(ctg) <- paste0(envrmt$path_lidar_norm,"{ID}","_norm_height")
  ctg <- lidR::normalize_height(ctg,lidR::knnidw())
  
  #- calculate chm
  lidR::opt_output_files(ctg) = paste0(envrmt$path_lidar_raster,"chm/{ID}","_chm")
  chm = rasterize_canopy(ctg, res = res, lidR::dsmtin())
  ker <- matrix(1,3,3)
  cchm <- raster::focal(chm, w = ker, fun = min, na.rm = TRUE)
  #chm_tif = get_vrt_img("grid_canopy",paste0(envrmt$path_lidar_raster,"chm/"),"chm")
  ttops_chm = locate_trees(cchm, lmf(8))
  algo1 = li2012(hmin = 5, R = 5)
  algo2 = lidR::dalponte2016(cchm,treetops = ttops_chm,th_tree = 5)
  algo3 = silva2016(chm, ttops_chm,exclusion = 0.5)
  
  
  #--- tree segmentation based on the point cloud using catalog
  # source: https://gis.stackexchange.com/questions/364546/individual-tree-segmentation-workflow-with-lidr-package
  #- tree segmentation using the 99 percentile filter
  opt_output_files(ctg) = paste0(envrmt$path_lidar_ctg,"hull_algo1/HULL_{XCENTER}_{YCENTER}")
  #opt_chunk_size(ctg) = 0
  opt_filter(ctg) <- "filter_noise(ctg, sensitivity = 1.2, res = 4)"
  ctg = segment_trees(ctg,  algo1 , uniqueness = "bitmerge")
  #laz = do.call(rbind, lapply(list.files(paste0(envrmt$path_lidar_raster,"hull_algo1/",prefix), pattern = "*.laz", full.names = TRUE),readLAS))
  saveRDS(ctg,ctgFN)
  
  #- tree filtering and calculate hull 
  opt_filter(ctg) <- "!is.na(treeID)"
  ctg@output_options$drivers$Spatial$extension = ".shp"
  opt_output_files(ctg) = paste0(envrmt$path_lidar_vector,"hull_sapflow/HULL_sapflow_{XCENTER}_{YCENTER}")
  hulls = catalog_apply(ctg=ctg, FUN = tree_fn)
  
  #- merge shapefiles
  #hulls <- list.files(paste0(envrmt$path_l_raster,"/",prefix), pattern = "*.shp", full.names = TRUE)
  seg = do.call(rbind, lapply(hulls, read_sf))
  st_write(seg,treesFN,append=F)
  #plot(st_geometry(seg))
  
  #- tmap plot
  tmap_mode("view")
  tm_shape(seg) + tm_fill(col = "zq95") 
  mapview(seg,zcol="zq95", fgb = FALSE)
  
  #- pixel_metrics 
  # NOTE: the metrics is performed on the manipulated (see above) ctg! 
  # So the results will differ from the raw lasfile/ctg
  sapflow_metrics <- pixel_metrics(ctg, .stdmetrics, res = res)
  raster::writeRaster(sapflow_metrics,sapflow_metricsFN)
  tmap_mode("view")
  mapview(sapflow_metrics) 
  
} else {
  dem = raster(demFN)
  chm = raster(chmFN)
  ctg = readRDS(ctgFN)
  trees = st_read(treesFN)
  sapflow_metrics = raster::stack(sapflow_metricsFN)
}