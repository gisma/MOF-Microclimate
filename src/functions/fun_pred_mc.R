  # functions for prediction

pred_rf_temp <- function(no){
  
  
  use_sp()
  const <- data.frame(
    #       doy = mx[no,2],
    #      hour = mx[no,3],
    Ta_10m = clim_stat$Ta_10m[no],
    rad_sw_in =clim_stat$rad_sw_in[no],
    Wind_direction_10m =clim_stat$Wind_direction_10m[no]  
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
                          const=const,progress="text")
  terra::writeRaster(pred,
                     paste0(envrmt$path_data_lev2,file_name,"_pred_rf_temp_",yr,".tif"),
                     overwrite=TRUE,gdal=c("COMPRESS=DEFLATE"))
  
}


pred_rf <- function(no){
  
  
  
  const <- data.frame(
    Huma_2m = clim_stat$Huma_2m[no]
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

