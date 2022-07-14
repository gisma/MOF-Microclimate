get_vrt_img = function(name,path,tag,raster = FALSE, stars= FALSE, terra = FALSE, fn = TRUE) {
  fList = list.files(path, 
                     pattern = glob2rx(paste0("*",tag,".tif")),
                     full.names = TRUE)
  
  res = invisible(gdalUtils::mosaic_rasters(gdalfile = fList,
                                            output_Raster = TRUE,
                                            dst_dataset = paste0(path,name,".tif"),verbose = FALSE))
  if (fn) return(paste0(path,name,".tif"))
  if (raster) return(raster::stack(paste0(path,name,".tif")))
  if (stars) return(read_stars(paste0(path,name,".tif")))
  if (terra) return(rast(paste0(path,name,".tif")))
  
}

template_raster = function(df,epsg){
  grid.MOF <- expand.grid(x = seq(from = st_bbox(bbox)$xmin,
                                  to = st_bbox(bbox)$xmax,
                                  by = 1),
                          y = seq(from = st_bbox(bbox)$ymin,
                                  to = st_bbox(bbox)$ymax,
                                  by = 1))
 
   # raster
  template_raster_MOF <-grid.MOF %>%
    raster::rasterFromXYZ(
      crs = epsg)
  return(template_raster_MOF)
}