# combine a list of whatever
laz = do.call(rbind, lapply(list.files(paste0(envrmt$path_level2,"/",prefix), pattern = "*.laz", full.names = TRUE),readLAS))
seg = do.call(rbind, lapply(hulls, sf::read_sf))
# list files 
hulls = list.files("path", pattern = "*.shp", full.names = TRUE)

# ---- voronoi tesselation of the tree stands ----
t_clust <- vect(tree_clust_sf)
# calculate voronoi
t_voronoi <- terra::voronoi(t_clust)
# reconvert it
sf_voronoi=st_as_sf(t_voronoi)
# mapview(sf_voronoi,zcol="cluster")
sf::st_write(sf_voronoi,voronoi_sf_FN, append= FALSE)
saveRDS(sf_voronoi,voronoi_rds_FN)

## ---- final voronoi tree clusters ----
# depreceated use tree clusters from envimetR/src/90_tree_cluser_analysis
# tree_clust_sf <- st_read(file.path(envrmt$path_level2,"sapflow_tree_all_cluster_sf.gpkg"))

tree_clust_rast <- gdalUtils::gdal_rasterize(src_datasource = voronoi_sf_FN,
                                             dst_filename =  voronoi_tif_FN,
                                             a = c("cluster"),
                                             l = "sf_voronoi",
                                             a_nodata = 0,
                                             a_srs = proj4,
                                             tr = c(res*1,res*1),
                                             te = c(chm@extent[1]  , 
                                                    chm@extent[3],
                                                    chm@extent[2],
                                                    chm@extent[4]),
                                             output_Raster = TRUE
)

tree_clust_rast <- rast(voronoi_tif_FN)
# plot(tree_clust_rast)

