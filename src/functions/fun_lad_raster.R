# ----   and more  ----
#fn = list.files(paste0(envrmt$path_lidar_ctg,"hull_algo3"),pattern = "HULL",full.names = TRUE)[[26]]

calc_specific_LAD_values = function(fList = NULL, raster=TRUE, visualize= FALSE){
  library(canopyLazR)
  newlad = newmetr = list()
  i=1
  for (fn in fList  ){
    
    laz.data <- laz.to.array(laz.file.path = fn, 
                             voxel.resolution = res, 
                             z.resolution = voxel.height,
                             use.classified.returns = FALSE)
    #saveRDS(laz.data,"all.rds")
    
    # Level the voxelized array to mimic a canopy height model
    level.canopy <- canopy.height.levelr(lidar.array = laz.data)
    cat("level canpoy: ",level.canopy$z.bin,"\n")
    if (length(level.canopy$z.bin) > htmin ){
      # Estimate LAD for each voxel in leveled array
      lad.estimates <- machorn.lad(leveld.lidar.array = level.canopy, 
                                   voxel.height = voxel.height, 
                                   beer.lambert.constant = 0.87)
      
      # Convert the LAD array into a single raster stack
      lad.raster <- lad.array.to.raster.stack(lad.array = lad.estimates, 
                                              laz.array = laz.data, 
                                              epsg.code = epsg_number)
      
      # Create a single LAI raster from the LAD raster stack
      lai.raster <- raster::calc(lad.raster, fun = sum, na.rm = TRUE)
      
      
      # Convert the list of LAZ arrays into a ground and canopy height raster
      grd.can.rasters <- array.to.ground.and.canopy.rasters(laz.data, epsg_number)
      
      # Calculate max LAD and height of max LAD
      max.lad <- lad.ht.max(lad.array = lad.estimates, 
                            laz.array = laz.data, 
                            ht.cut = htmin, 
                            epsg.code = epsg_number)
      
      # Calculate the ratio of filled and empty voxels in a given column of the canopy
      empty.filled.ratio <- canopy.porosity.filled.ratio(lad.array = lad.estimates,
                                                         laz.array = laz.data,
                                                         ht.cut = htmin,
                                                         epsg.code = epsg_number)
      
      # Calculate the volume of filled and empty voxles in a given column of the canopy
      empty.filled.volume <- canopy.porosity.filled.volume(lad.array = lad.estimates,
                                                           laz.array = laz.data,
                                                           ht.cut = htmin,
                                                           xy.res = res,
                                                           z.res = voxel.height,
                                                           epsg.code = epsg_number)
      
      # Calculate the within canopy rugosity
      within.can.rugosity <- rugosity.within.canopy(lad.array = lad.estimates,
                                                    laz.array = laz.data,
                                                    ht.cut = htmin,
                                                    epsg.code = epsg_number)
      
      # Calculate the heights of various LAD quantiles
      ht.quantiles <- lad.quantiles(lad.array = lad.estimates,
                                    laz.array = laz.data,
                                    ht.cut = htmin,
                                    epsg.code = epsg_number)
      
      # Calculate various canopy volume metrics from Lefsky
      can.volume <- canopy.volume(lad.array = lad.estimates,
                                  laz.array = laz.data,
                                  ht.cut = htmin,
                                  xy.res = res,
                                  z.res = voxel.height,
                                  epsg.code =epsg_number)
      
      # We can calculate the depth of the euphotic zone by dividing by the volume of the voxel
      euphotic.depth <- can.volume$euphotic.volume.column.raster / ( res * res *voxel.height)
      
      # Calculate the top of canopy rugosity volume
      toc.rugos <- toc.rugosity(chm.raster = grd.can.rasters$chm.raster,
                                xy.res = res,
                                z.res = voxel.height)
      
      # Plot the lai raster
      #plot(lai.raster)
      # 
      # # Plot the ground raster
      # plot(grd.can.rasters$ground.raster)
      # 
      # # Plot the canopy height raster
      # plot(grd.can.rasters$canopy.raster)
      # 
      # # Plot the canopy height model raster
      # plot(grd.can.rasters$chm.raster)
      # 
      # # Plot the max LAD raster
      # plot(max.lad$max.lad.raster)
      # 
      # # Plot the height of max LAD raster
      # plot(max.lad$max.lad.ht.raster)
      # 
      # # Plot filled voxel ratio raster
      # plot(empty.filled.ratio$filled.raster)
      # 
      # # Plot porosity voxel ratio raster
      # plot(empty.filled.ratio$porosity.raster)
      # 
      # # Plot filled voxel volume raster
      # plot(empty.filled.volume$filled.raster)
      # 
      # # Plot porosity voxel volume raster
      # plot(empty.filled.volume$porosity.raster)
      # 
      # # Plot the standard deviation of LAD within a vertical column raster
      # plot(within.can.rugosity$vertical.sd.lad.raster)
      # 
      # # Plot within canopy rugosity
      # plot(within.can.rugosity$rugosity.raster)
      # 
      # # Plot the height of the 10th quantile
      # plot(ht.quantiles$quantile.10.raster)
      # 
      # # Plot the height of the 25th quantile
      # plot(ht.quantiles$quantile.25.raster)
      # 
      # # Plot the height of the 50th quantile
      # plot(ht.quantiles$quantile.50.raster)
      # 
      # # Plot the height of the 75th quantile
      # plot(ht.quantiles$quantile.75.raster)
      # 
      # # Plot the height of the 90th quantile
      # plot(ht.quantiles$quantile.90.raster)
      # 
      # # Plot the height of the mean LAD
      # plot(ht.quantiles$mean.raster)
      # 
      # # Plot the volume of the euphotic zone for each column
      # plot(can.volume$euphotic.volume.column.raster)
      # 
      # # Plot the total leaf area in the euphotic zone for each column
      # plot(can.volume$euphotic.tla.column.raster)
      # 
      # # Plot the depth of the euphotic zone
      # plot(euphotic.depth)
      # 
      # # Plot the volume of the oligophotic zone for each column
      # plot(can.volume$oligophotic.volume.column.raster)
      # 
      # # Plot the total leaf area in the oligophotic zone for each column
      # plot(can.volume$oligophotic.tla.column.raster)
      # 
      # # Plot the volume of the empty space within a given colume
      # plot(can.volume$empty.volume.column.raster)
      # 
      # # Plot the volume of the empty space within a 3x3 moving window
      # plot(can.volume$empty.canopy.volume.raster)
      # 
      # # Plot the volume of the euphotic zone within a 3x3 moving window
      # plot(can.volume$filled.canopy.euphotic.raster)
      # 
      # # Plot the volume of the oligophotic zone within a 3x3 moving window
      # plot(can.volume$filled.canopy.oligophotic.raster)
      # 
      # # Plot the total leaf area of the euphotic zone within a 3x3 moving window
      # plot(can.volume$filled.canopy.euphotic.tla.raster)
      # 
      # # Plot the total leaf area of the oligophotic zone within a 3x3 moving window
      # plot(can.volume$filled.canopy.oligophotic.tla.raster)
      # 
      # # Plot the top of canopy rugosity volume
      # plot(toc.rugos)
      lad_names = c("toc.rugos","filled.canopy.oligophotic.tla","filled.canopy.euphotic.tla",
                    "filled.canopy.oligophotic","filled.canopy.euphotic","empty.canopy.volume",
                    "empty.volume.column","oligophotic.tla.column","oligophotic.volume.column",
                    "euphotic.depth","euphotic.tla.column","euphotic.volume.column",
                    "quantiles.mean","quantile.90","quantile.75","quantile.50","quantile.25","quantile.10",
                    "within.can.rugosity","within.can.rugosity.vertical.sd.lad","empty.filled.volume.porosity",
                    "empty.filled.volume.filled","empty.filled.ratio.porosity","empty.filled.ratio.filled",
                    "max.lad.max.lad.ht","max.lad.max.lad","grd.can.rasters.chm","grd.can.rasters.canopy",
                    "grd.can.rasters.ground","lai.raster")
      newmetr[[i]] =  raster::stack(toc.rugos,can.volume$filled.canopy.oligophotic.tla.raster,
                                    can.volume$filled.canopy.euphotic.tla.raster,
                                    can.volume$filled.canopy.oligophotic.raster,
                                    can.volume$filled.canopy.euphotic.raster,
                                    can.volume$empty.canopy.volume.raster,
                                    can.volume$empty.volume.column.raster,
                                    can.volume$oligophotic.tla.column.raster,
                                    can.volume$oligophotic.volume.column.raster,
                                    euphotic.depth,
                                    can.volume$euphotic.tla.column.raster,
                                    can.volume$euphotic.volume.column.raster,
                                    ht.quantiles$mean.raster,
                                    ht.quantiles$quantile.90.raster,
                                    ht.quantiles$quantile.75.raster,
                                    ht.quantiles$quantile.50.raster,
                                    ht.quantiles$quantile.25.raster,
                                    ht.quantiles$quantile.10.raster,
                                    within.can.rugosity$rugosity.raster,
                                    within.can.rugosity$vertical.sd.lad.raster,
                                    empty.filled.volume$porosity.raster,
                                    empty.filled.volume$filled.raster,
                                    empty.filled.ratio$porosity.raster,
                                    empty.filled.ratio$filled.raster,
                                    max.lad$max.lad.ht.raster,
                                    max.lad$max.lad.raster,
                                    grd.can.rasters$chm.raster,
                                    grd.can.rasters$canopy.raster,
                                    grd.can.rasters$ground.raster,
                                    lai.raster)
      names(newmetr[[i]]) = lad_names
      newlad[[i]] = lad.raster  
      plot(newmetr[[i]])
    }
    i = i + 1
  }
  
  if (raster){
    rasters1.mosaicargs <- newmetr
    rasters1.mosaicargs$fun <- mean
    m <- do.call(mosaic, rasters1.mosaicargs)  
    names(m) = lad_names
    return(m)
  }
  if (visible){
    names(m) = lad_names
    return(mapvieww::mapview(m))
  }
  
  
}