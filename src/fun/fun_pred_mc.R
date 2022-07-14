# function for prediction

pred_mc <- function(doy,
                    hour,
                    model){
  hour_chr <-  sprintf("%02d", hour)
  
  file_name <- paste("tot_rad_", doy, "_", hour_chr, ".00", sep = "")
  
  raster <- raster(readRAST(file_name))
  
  raster <- crop(raster, pred_stack)
  
  index = which(preProcValues$method$center == "rad")
  
  raster = (raster - preProcValues$mean[index])/preProcValues$std[index]

  pred_stack2 <- stack(pred_stack_scaled, raster)
  names(pred_stack2) <- c("dtm", "slope", "aspect", "TPI","CHM", "ip_ground", "cluster","rad")
  pred_stack2 = dropLayer(pred_stack2, 7)
  
  hour_scaled = as.numeric((hour - preProcValues$mean[which(preProcValues$method$center == "hour")])/preProcValues$std[which(preProcValues$method$center == "hour")])
  doy_scaled = as.numeric((doy - preProcValues$mean[which(preProcValues$method$center == "doy")])/preProcValues$std[which(preProcValues$method$center == "doy")])
  
  i <- which(clim_stat$doy == doy & clim_stat$hour == hour)
  j <- which(trainDF$doy == doy_scaled & trainDF$hour == hour_scaled)
  j <- j[1]
  
  const <- data.frame(doy = doy, 
                      hour = hour,
                      Ta_10m = clim_stat$Ta_10m[i],
                      #Windspeed_2m = clim_stat$Windspeed_2m[i],
                      Windspeed_2m = clim_stat$Windspeed_2m[i],
                      Wind_direction_10m = clim_stat$Wind_direction_10m[i],
                      Wind_direction_3m = clim_stat$Wind_direction_3m[i],
                      rad_sw_in = clim_stat$rad_sw_in[i],
                      rad_lw_in = clim_stat$rad_lw_in[i],
                      #Ta_2m = clim_stat$Ta_2m[i],
                      Huma_10m = clim_stat$Huma_10m[i],
                      cloudiness = trainDFclean$cloudiness[j],
                      rad_Klimastation = trainDFclean$rad_Klimastation[j],
                      day = trainDFclean$day[j])
  
  index = vector()
  
  for (i in 1:length(names(const))) {
    if (names(const)[i] %in% preProcValues$method$center) {
      index[i] = which(preProcValues$method$center == names(const)[i])
    } else {
      index[i] = NA
    }
  }
  
  const_scaled = (const - preProcValues$mean[index])/preProcValues$std[index]
  
  #f = list(levels(trainDF$cluster))
  pred_scaled <- predict(pred_stack2, model, const = const_scaled, allow.new.levels = TRUE, progress = "text")
  
  pred = (pred_scaled * preProcValues$std[which(preProcValues$method$center == "temp")]) + preProcValues$mean[which(preProcValues$method$center == "temp")]

  }
