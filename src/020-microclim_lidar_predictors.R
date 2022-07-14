#
# Name: 20_microclimate_predictor_stack.R
# Type: control script 
# Description:  creates a predictor stack derive raster/tree segment based indices
#               for microclimate modelling
#               canopy height model, digital elevation models and standard statistics
#         - irridiance (r.sun.hourly)
#         - topographic information (based on dem): dem, dsm , slope, aspect, TPI
#         - forest structure (based on lidar): chm
#                            tree cluster based LAD and std-metrics
#                            gridmetrics including VCI LAD entropy, ipground std-metrics 
# Input: las_file 
# Output: - las catalog 
#         - raster grids and vector data sets as mentioned
#         - prediction stack with all layers 
#               forest_structure_metrics
#               trees_structure_metrics 
#               trees_lad_DF 
# Copyright: GPL (>= 3),Lena Perzlmaier, Chris Reudenbach, creuden@gmail.com 
# Date: 2022-01-24
#

# ----project setup----
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# Load the library
appendpackagesToLoad = c("lidR","future","lwgeom","tmap","canopyLazR")
appendProjectDirList =  c("data/lidar/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/pmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/tmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/topo/")
source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)

# ---additional variables----

# filenames 
ctg_baseFN= paste0(envrmt$path_lidar_ctg,"base/")
treesFN = paste0(envrmt$path_lidar_vector,"LIR0hmin5_segmentation.gpkg")
las_fileFN = paste0(envrmt$path_lidar_orig,"MOF_lidar_2018.las")
topoFN = paste0(envrmt$path_topo,"pred_topo.tif")
pmetricsFN = paste0(envrmt$path_lidar_raster,"pmetrics/all_pixel_metrics.tif")
pLADFN = paste0(envrmt$path_lidar_raster,"pmetrics/pixel_LAD_metrics.tif")
pentropyFN = paste0(envrmt$path_lidar_raster,"pmetrics/pixel_entropy_metrics.tif")
pVCIFN = paste0(envrmt$path_lidar_raster,"pmetrics/pixel_VCI_metrics.tif")
ladFN = paste0(envrmt$path_lidar_vector,"tmetrics/lad.gpkg")
ladhullFN = paste0(envrmt$path_vector,"lad_hull.gpkg")
ladhull_metricsFN = paste0(envrmt$path_tmetrics,"lad_hull_raster.tif")
treelad_cleanFN = paste0(envrmt$path_lidar,"trees_lad_clean.rds")
best_clu_no_FN = paste0(envrmt$path_lidar,"opticlust.rds")
tree_clus_sfFN = file.path(paste0(envrmt$path_tmetrics,"tree_cluster.gpkg"))
tree_clus_rasFN = file.path(paste0(envrmt$path_tmetrics,"tree_cluster.tif"))
predstack_forest_metricsFN = file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif")

# calculate or just read final data set
calculate = TRUE

# in addition (re-)calculate tree hull based analysis
tree_hulls = TRUE

# in addition (re-)calculate GRASS irridiance calculation
sol = FALSE

#- resolution for catalog operations in meter
res = 1.0
voxel.height = 1.0
zmax = 40
htmin = 4.0
k= 0.87

#- parallelisation number of Cores
future::plan(multisession, workers = 24)
set_lidr_threads(24)

# ---- start processing ----

if (calculate){
  # ---- define first basic lidR catalog ----
  ctg <- lidR::readLAScatalog(las_fileFN)
  projection(ctg) <- epsg_number
  lidR::opt_chunk_size(ctg) = 320
  lidR::opt_chunk_buffer(ctg) <- 20
  lidR::opt_progress(ctg) <- TRUE
  lidR::opt_laz_compression(ctg) <- TRUE
  lidR::opt_filter(ctg) = "filter_noise(ctg, sensitivity = 1.2, res = res)"
  ctg@output_options$drivers$Raster$param$overwrite <- TRUE
  ctg@output_options$drivers$SpatRaster$param$overwrite <- TRUE
  ctg@output_options$drivers$SpatVector$param$overwrite <- TRUE
  ctg@output_options$drivers$stars$param$overwrite <- TRUE
  ctg@output_options$drivers$sf$param$overwrite <- TRUE
  ctg@output_options$drivers$Spatial$param$overwrite <- TRUE
  ctg@output_options$drivers$LAS$param$overwrite <- NULL
  lidR::opt_output_files(ctg) = paste0(envrmt$path_lidar_ctg,"base/{ID}","_base")
  ctg_base = lidR::normalize_height(ctg,lidR::knnidw(k = 6L, p = 2))
  
  # ---- dsm, dem chm and basic morphometrics ----
  # source: https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
  lidR::opt_output_files(ctg) <- paste0(envrmt$path_lidar_raster,"dem/{ID}","_dem")
  dem <- rasterize_terrain(ctg, res = res,lidR::knnidw(k = 6L, p = 2))
  
  #- dsm 
  lidR::opt_output_files(ctg) <- paste0(envrmt$path_lidar_raster,"dsm/{ID}","_dsm") # add output filname template
  dsm = lidR::rasterize_canopy(ctg, res = res, algorithm = pitfree())
  
  #- chm 
  # NOTE chm is calculated on ctg_base
  lidR::opt_output_files(ctg_base) = paste0(envrmt$path_lidar_raster,"chm/{ID}","_chm")
  chm = rasterize_canopy(ctg_base, res = res,pitfree(c(0,2,5,10,15)) ,pkg = "terra")
  
  #-  standard morphometric parameters based on dsm
  # slope
  slope <- terrain(dem, "slope")
  # aspect 
  aspect <- terrain(dem,  "aspect", unit = "degrees")
  # topographic position index 
  TPI <- terrain(dsm, "TPI")
  ker <- matrix(1,3,3)
  TPI <- terra::focal(TPI, w = ker, fun = mean, na.rm = TRUE)
  TPI[TPI < 0 ] <- - 1
  TPI[TPI > 0 ] <-  1
  
  # stack basic topographic information
  topo <- c(dem, dsm, chm, slope, aspect, TPI)
  names(topo) <- c("dem","dsm","chm", "slope", "aspect", "TPI")
  plot(topo)
  terra::writeRaster(topo,  topoFN , overwrite = TRUE)
  rm(dem,dsm,chm, slope, aspect, TPI)
  
  # ---- pixel_metrics----
  # stdmetrics
  ctg_base <- lidR::readLAScatalog(ctg_baseFN)
  ctg_base@output_options$drivers$SpatRaster$param$overwrite <- TRUE
  ctg_base@output_options$drivers$sf$param$overwrite <- TRUE
  opt_output_files(ctg_base) = paste0(envrmt$path_pmetrics,"std/{ID}","_base_metrics")
  pixel_stdmetrics = pixel_metrics(ctg_base, func = .stdmetrics, res = res)
  metrics_FN = get_vrt_img("pixel_metrics_MOF",paste0(envrmt$path_pmetrics,"std"),"base_metrics")
  pixel_stdmetrics =rast(metrics_FN)
  names(pixel_stdmetrics) = stm_names[1:56]  
  plot(pixel_stdmetrics)
  
  
  # cov_LAD accu 
  opt_output_files(ctg_base) = paste0(envrmt$path_pmetrics,"lad/{ID}","_LAD_metrics")
  pixel_LAD = pixel_metrics(ctg_base, ~as.numeric(cv(LAD(Z, dz = 1, k= k)$lad)), res)
  lad_FN = get_vrt_img("LAD_metrics_MOF",paste0(envrmt$path_pmetrics,"lad/"),"LAD_metrics")
  pixel_LAD =rast(lad_FN)
  plot(pixel_LAD)
  
  
  # entropy accu
  opt_output_files(ctg_base) = paste0(envrmt$path_pmetrics,"entropy/{ID}","_entropy_metrics")
  pixel_entropy = pixel_metrics(ctg_base, ~as.numeric(entropy(Z, by = voxel.height)), res)
  entropy_FN = get_vrt_img("entropy_metrics_MOF",paste0(envrmt$path_pmetrics,"entropy/"),"entropy_metrics")
  pixel_entropy = rast(entropy_FN)
  plot(pixel_entropy)
  
  # VCI accu
  opt_output_files(ctg_base) = paste0(envrmt$path_pmetrics,"vci/{ID}","_VCI_metrics")
  pixel_VCI = pixel_metrics(ctg_base, ~as.numeric(VCI(Z, zmax = zmax, by = voxel.height)), res)
  VCI_FN = get_vrt_img("VCI_metrics_MOF",paste0(envrmt$path_pmetrics,"vci/"),"VCI_metrics")
  pixel_VCI = rast(VCI_FN)
  plot(pixel_VCI)
  
  # ipground  percentage of intensity returned by points classified as "ground" class == 2
  # extracted from lidR sourcecode  ipground <- c(metrics, list(ipground = sum(i[class == 2])/itot*100))
  lidR::opt_output_files(ctg_base) = paste0(envrmt$path_pmetrics,"ip_ground/{ID}","_ip_ground")  
  ip_ground <- pixel_metrics(ctg_base, ~sum(Intensity[Classification == 2])/sum(Intensity)*100, res = res)
  ipgroundFN = get_vrt_img("ip_ground_MOF",paste0(envrmt$path_pmetrics,"ip_ground/"),"ip_ground")
  pixel_ip_ground = rast(ipgroundFN)
  plot(ip_ground)
  
  # ---- stack pixelbased metrics ----
  pmetrics = c(pixel_VCI,pixel_entropy,pixel_LAD,pixel_ip_ground,pixel_stdmetrics)
  terra::writeRaster(pmetrics,pmetricsFN,overwrite=TRUE)
  pmetrics_all = rast(pmetricsFN)
  plot(pmetrics_all)
  
  rm(pmetrics,pixel_VCI,pixel_entropy,pixel_LAD,pixel_stdmetrics,pixel_ip_ground)
  
  if (tree_hulls){
    
    # ---- tree based analysis  ----
    # segmentation using Li et al. 
    # for chm based approach see  https://gis.stackexchange.com/questions/364546/individual-tree-segmentation-workflow-with-lidr-package
    ctg_base <- lidR::readLAScatalog(ctg_baseFN)
    ctg_base@output_options$drivers$SpatRaster$param$overwrite <- TRUE
    ctg_base@output_options$drivers$sf$param$overwrite <- TRUE
    opt_stop_early(ctg_base) <- FALSE
    opt_output_files(ctg_base) = paste0(envrmt$path_lidar_ctg,"segmented/seg_{ID}")
    ctg_seg = segment_trees(ctg_base,  li2012() , uniqueness = "bitmerge")
    
    
    # ---- calculate hulls ----
    # NOTE hulls are calculated on ctg_seg
    opt_filter(ctg_seg) = "!is.na(treeID)"
    ctg_seg@output_options$drivers$Spatial$extension = ".gpkg"
    opt_output_files(ctg_seg) = paste0(envrmt$path_lidar_vector,"tree_segs/LIR0hmin5_hulls")
    hulls = catalog_apply(ctg_seg, FUN = tree_fn) # tree_fn is a user defined function
    # - combine tree segements getting the LADs----
    seg = do.call(rbind, lapply(hulls, sf::read_sf))
    st_write(seg,treesFN,append=FALSE)
    tree_hulls = st_read(treesFN)
    
    # call adapted function from leafR to calculate LAD Profiles
    # http://dx.doi.org/10.3390/rs11010092   
    LAD_voxDF =  lad.voxels(ctg_base,grain.size = 1, k=k,maxP = zmax)
    lad_sf = st_as_sf(LAD_voxDF,coords = c("X", "Y"), crs = epsg_number)
    sf::st_write(lad_sf,ladFN)
    lad_sf = st_read(ladFN)
    # join LAD & hulls
    stat_hulls = st_intersection(lad_sf,tree_hulls)
    
    stat_hulls = st_write(stat_hulls,ladhullFN,append = FALSE)
    stat_hulls = st_read(ladhullFN)
    # create empty raster
    stdmetrics_r = template_raster(bbox,epsg_number)
    # rasterize the data
    stdmetrics_r = fasterize::fasterize(st_buffer(stat_hulls,dist = 1.),
                                        field = names(stat_hulls), stdmetrics_r,fun="sum")
    #stdmetrics_r = raster::rasterize(as(stat_hulls,"Spatial"), stdmetrics_r, fun = mean)
    terra::writeRaster(stdmetrics_r,ladhull_metricsFN,overwrite=TRUE)
    trees_lad_stdmetrics = rast(ladhull_metricsFN)
    rm(stdmetrics_r)
    # drop geometry
    trees_lad = st_drop_geometry(stat_hulls)
    tree_ID = stat_hulls$treeID
    trees_lad = trees_lad[keep_names]
    
    ##--- remove rows with NA
    trees_lad = trees_lad[complete.cases(trees_lad) ,]
    
    # 
    # drop empty cols
    trees_lad=trees_lad[ , colSums(is.na(trees_lad)) == 0]
    trees_lad=trees_lad[ , colSums(trees_lad) > 0]
    ##--- filter correlations that are > cor_cutoff
    filt = findCorrelation(cor(trees_lad), cutoff = 0.9)
    trees_lad = trees_lad[,-filt]
    
    # save the clean data frame containing the lidar derived statistics
    saveRDS(trees_lad,treelad_cleanFN)
    trees_lad = readRDS(treelad_cleanFN)
    
    
    # ---- cluster analysis ----
    #- prepare data table
    # create header to drop  fantasy tree_hulls above 44.5 meters
    #header_height=paste0("z_",seq(44.5,100.5,1))
    
    # define df for clustering without the following attributes
    data_clust = trees_lad[ ,!(names(trees_lad) %in% c("treeID","x","y","n","area"))]
    
    # perform PCA for much faster clustering with identical results
    pca_dat = stats::princomp(data_clust)$scores[, 1:10]
    # brute force data sampling for cluster number estimate
    data = sample_n(data_clust,floor(nrow(data_clust)* 0.01))
    
    # highspeed analysis of how many clusters are meaningful
    # all clusternumbers below threshold 0.85 are ok but prefer stable results
    opt = ClusterR::Optimal_Clusters_KMeans(data, max_clusters = 25, plot_clusters = T,num_init = 50,
                                            criterion = 'distortion_fK',
                                            initializer = 'kmeans++',
                                            seed = 123)
    saveRDS(opt,best_clu_no_FN)
    #- clustering
    # kmeans clustering with the number of clusters as derived by opt
    # this is difficult to automize due to the effect that I prefer the cutoff level
    # after the first stable grow of the threshold value the basically the best value 
    # is two due to build up and non build up areas than we see in all runs a stable 
    # grow up to roughly < 15 clusters with a values > .85  clusters
    clust_model = ClusterR::KMeans_arma(data_clust, clusters = 10, n_iter = 500,
                                        seed_mode = "random_subset",
                                        verbose = T, CENTROIDS = NULL)
    
    #- prediction on all data
    trees_lad$cluster = as.integer(ClusterR::predict_KMeans(data_clust, clust_model))
    #trees_lad$cluster= as.double(trees_lad$cluster)
    
    # join the results to sf hulls polygons
    t_cluster = inner_join(tree_hulls,trees_lad[,names(trees_lad) %in% keep_lad_names,])
    
    #t_cluster = t_cluster[,  !(names(t_cluster) %in% "geometry")]
    #t_cluster = t_cluster[complete.cases(t_cluster) ,]
    #tree_clust_sf = st_as_sf(t_cluster,coords = c("x", "y"), crs = epsg_number)
    sf::st_write(t_cluster,tree_clus_sfFN, append= FALSE)
    
    
    tree_clust_rast <- gdalUtils::gdal_rasterize(src_datasource = tree_clus_sfFN,
                                                 dst_filename = tree_clus_rasFN,
                                                 a = c("cluster"),
                                                 a_nodata = 0,
                                                 a_srs = proj4,
                                                 tr = c(res,res),
                                                 te = c(st_bbox(bbox)$xmin,  st_bbox(bbox)$ymin,
                                                        st_bbox(bbox)$xmax,  st_bbox(bbox)$ymax))
    tree_clust_ras <- rast(tree_clus_rasFN)
    plot(as(tree_clust_ras,"Raster"))
    # tc=c(tree_clust_ras,trees_lad_stdmetrics)
    # terra::writeRaster(tc)
    # 
  }
  # ---- stack forest structure metrics----
  pindices = rast(topoFN)
  pmetrics <- rast(pmetricsFN)
  tree_clust_ras <- rast(tree_clus_rasFN)
  
  forest_structure_metrics <-c(pindices,tree_clust_ras,pmetrics)
  #names(forest_structure_metrics)[12:length(names(forest_structure_metrics))] <- stm_names[1:56]
  # save forest structure stack
  terra::writeRaster(forest_structure_metrics, predstack_forest_metricsFN,overwrite=T)
  
  ## ---- calculate irradiance----
  if (sol)  {
    # NOTE: the GRASS based calxculations are utilized via a wrapper logic
    # the resulting files are NOT directly available through GDAL driver logic 
    link2GI::linkGRASS7(dem, gisdbase = root_folder, location = "MOF2")
    
    # import dtm
    rgrass7::execGRASS("r.in.gdal", 
                       flags = "o",
                       parameters = list(
                         input = paste0(envrmt$path_lidar_raster,"dem/rasterize_terrain.vrt"),
                         output = "dtm"
                       ))
    
    rgrass7::execGRASS("r.info", 
                       parameters = list(map = "dtm"))
    
    # import topo pred
    rgrass7::execGRASS("r.in.gdal", 
                       flags = "o",
                       parameters = list(
                         input = file.path(envrmt$path_topo,"/pred_topo.tif"),
                         output = "pred_topo"
                       ))
    
    rgrass7::execGRASS("r.info", 
                       parameters = list(map = "pred_topo.1"))
    
    # set extent of the region to the extent of the dtm
    rgrass7::execGRASS("g.region", 
                       parameters = list(raster = "dtm"))
    
    # install r.sun.hourly
    rgrass7::execGRASS("g.extension",
                       parameters = list(extension = "r.sun.hourly"))
    
    
    # calculate hourly irradiance for the whole study period
    for (i in 1:1458) {
      
      output <- paste("tot_rad_", i, sep = "")
      rgrass7::execGRASS("r.sun.hourly", flags=c("overwrite", "quiet"), elevation="pred_topo.1", aspect="pred_topo.3", slope="pred_topo.2", mode="mode1",
                         start_time=0, end_time=23, day=i, year=2018, glob_rad_basename=output)
    }
    
  }   
} else {
  forest_structure_metrics <- rast(predstack_forest_metricsFN)  
  trees_structure_metrics = rast(ladhull_metricsFN)
  trees_lad_DF = readRDS(treelad_cleanFN)
}  
