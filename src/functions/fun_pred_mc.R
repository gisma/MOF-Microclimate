# function for prediction

pred_mc <- function(no){
  
  
  use_sp()
  const <- data.frame(
     doy = mx[no,2],
    hour = mx[no,3],
    Ta_10m = clim_stat$Ta_10m[no]
  )
  
  hour_chr <-  sprintf("%02d", mx[no,3])
  file_name <- paste("tot_rad_", mx[no,2], "_", hour_chr, ".00", sep = "")
  yr = lubridate::year(as.POSIXct(mx[no,1]))

  r <- raster::raster(rgrass7::read_RAST(file_name))
  tmp_stack <- raster::stack(pred_stack, r)
  names(tmp_stack) <- c(names(pred_stack),"rad")
  pred <- raster::predict(object = tmp_stack, 
                          model = model, 
                          allow.new.levels=TRUE, 
                          const=const)
  terra::writeRaster(pred,
                     paste0(envrmt$path_data_lev2,file_name,"_pred_inter_",yr,".tif"),
                     overwrite=TRUE,gdal=c("COMPRESS=DEFLATE"))
  
}


pred_rf <- function(no){
  
  

  const <- data.frame(
    doy = mx[no,2],
    hour = mx[no,3],
    Ta_10m = clim_stat$Ta_10m[no],
    Windspeed_10m = clim_stat$Windspeed_10m[no]
  )
  
  hour_chr <-  sprintf("%02d", mx[no,3])
  file_name <- paste("tot_rad_", mx[no,2], "_", hour_chr, ".00", sep = "")
  yr = lubridate::year(as.POSIXct(mx[no,1]))

  pred <- raster::predict(object = pred_stack, 
                          model = model_rf , 
                          allow.new.levels=TRUE, 
                          const=const)
  terra::writeRaster(pred,
                     paste0(envrmt$path_data_lev2,file_name,"_pred_rf_hum_",yr,".tif"),
                     overwrite=TRUE,gdal=c("COMPRESS=DEFLATE"))
  
}
