#------------------------------------------------------------------------------
# Name: make_CHM_Tiles.R
# Type: control script 
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  script creates a canopy height model from generic Lidar 
#              las data using the lidR package
# Data: regular las LiDAR data sets 
# Copyright:GPL (>= 3) 
# Date: 2021-12-10
#------------------------------------------------------------------------------
# 0 - specific setup
#-----------------------------
#-----------------------------
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()
# OPTIONAL additional packages
appendpackagesToLoad = c("lidR","future","lwgeom","tmap")
# OPTIONAL additinal folders
appendProjectDirList =  c("data/lidar/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/")

source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)

#-- Further customization of the setup by the user this section 
#-- can be freely customized only the definition of additional packages 
#-- and directory paths MUST be done using the two variables 
#-- appendpackagesToLoad and appendProjectDirList
#-- feel free to remove this lines if you do not need them
# define  additional packages comment if not needed

# MANDANTORY: calling the setup script also DO NOT change this line
source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)

# 1 - start script
#-----------------------------

# if you do not have original data you may uncomment the following lines and download a sample
#---- NOTE file size is about 12MB
# utils::download.file(url="https://github.com/gisma/gismaData/raw/master/uavRst/data/lidR_data.zip",
#                      destfile=paste0(envrmt$path_l_raw,"/chm.zip"))
# unzip(paste0(envrmt$path_l_raw,"/chm.zip"),
#       exdir = envrmt$path_l_raw,  
#       overwrite = TRUE)

#---- Get all *.las files of a folder into a list
las_files = list.files(envrmt$path_l_raw,
                       pattern = glob2rx("*.las"),
                       full.names = TRUE)

#--- create CHM with lidR catalog
#- source: https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
future::plan(multisession, workers = 6)
set_lidr_threads(6)
ctg <- lidR::readLAScatalog(las_files[[2]])
lidR::projection(ctg) <- 25832
lidR::opt_chunk_size(ctg) = 500
lidR::opt_chunk_buffer(ctg) <- 5
lidR::opt_progress(ctg) <- TRUE
lidR::opt_laz_compression(ctg) <- TRUE
ctg@output_options$drivers$Raster$param$overwrite <- TRUE

#--- height normalisation within the point cloud

lidR::opt_output_files(ctg) <- paste0(envrmt$path_l_norm,"/{ID}","_norm_height")
norm_ctg <- lidR::normalize_height(ctg,lidR::tin())

# reassign the projection
sp::proj4string(ctg) <- sp::CRS(proj4)

# calculate the chm with the pitfree algorithm
chm = lidR::rasterize_canopy(ctg, 0.25, pitfree(c(0,2,5,10,15), c(0,1), subcircle = 0.2))

# write it to tif
raster::writeRaster(chm,file.path(envrmt$path_l_raster,"mof_chm_one_tile.tif"),overwrite=TRUE) 


 # - visualize 
-------------------
  
# call mapview with some additional arguments
  mapview::mapview(chm,
          legend=TRUE, 
          layer.name = "canopy height model",
          col = mvTop(256),
          alpha.regions = 0.7)