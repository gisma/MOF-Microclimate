# --------------------------------------------------
# Name: 015-microclim_prepare_climate_station_data
# Author: Lena Perzlmaier Chris Reudenbach
# Description: Clean and prepare Grubenwiese climate station data 
#         This data serves as input for the model
#         The idea is to have a model, which needs only the measured 
#         temperature/relative humidity/cloudiness and date/time
#         in order to have a microclimate modell of the MOF as output
#         data will be filled and checked using hmisc
#         - period 01.01.2018 - 31.12.2021
#         - aggregated to hourly data
# Data:   downloaded from http://lcrs.geographie.uni-marburg.de/lcrs/data_pre.do?citid=302
# Output: "klimastation_wiese_hourly.rds"
# -------------------------------------------------


# 0 - load packages
#------------------
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()
appendpackagesToLoad = c("lidR","future","lwgeom","tmap","plyr","dplyr","tidyr","tidyverse","lubridate","oce","forecast","ggplot2","VIM","Hmisc","mice")
appendProjectDirList =  c("data/lidar/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/",
                          "data/lidar/MOF_lidar_2018/pmetrics/",
                          "data/lidar/MOF_lidar_2018/tmetrics/",
                          "data/lidar/MOF_lidar_2018/topo/",
                          "data/measurements/measure_org",
                          "data/measurements/measure_lev0",
                          "data/measurements/measure_lev1",
                          "data/measurements/")

source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)

gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")
gw_clim_csvFN = file.path(envrmt$path_measure_org, "302_KlimaturmCaldernWiese_complete.csv")

startdate = "2018-01-01"
enddate   = "2022-01-03"

# 2- start code
#------------------

clim_stat <- read.table(gw_clim_csvFN,
                        header = TRUE,
                        sep = ",",
                        dec = ".",
                        na.strings = "NULL")

# convert to posixct and filter to relevant time frame
clim_stat$datetime <-  as.POSIXct(clim_stat$datetime, format = "%Y-%m-%d %H:%M:%S")
clim_stat <- clim_stat %>% filter(between(datetime, as.POSIXct(startdate), as.POSIXct(enddate)))

plot(clim_stat$datetime, clim_stat$Ta_2m, type = "l")

# summarize data to hourly data
clim_stat$date_time_hourly <- substr(clim_stat$datetime, 1, 13)
clim_stat$date_time_hourly <- as.POSIXct(clim_stat$date_time_hourly, format = "%Y-%m-%d %H")
ts = seq.POSIXt(min(clim_stat$date_time_hourly,na.rm = TRUE),max(clim_stat$date_time_hourly,na.rm = TRUE),by="hour")
df=data.frame(date_time_hourly=ts)
clim_stat <- full_join(df,clim_stat)
clim_stat_2 = aggregate(clim_stat,  list(clim_stat$date_time_hourly),  mean)
 
plot(clim_stat_2$date_time_hourly, clim_stat_2$Ta_2m ,type = "l")

# control on gaps
md.pattern(clim_stat_2)
mice_plot <- aggr(clim_stat_2, col=c('navyblue','yellow'),
                    numbers=TRUE, sortVars=TRUE,
                    labels=names(clim_stat_2), cex.axis=.3,
                    gap=3, ylab=c("Missing data","Pattern"))
# fill gaps   
f <- aregImpute(~Ta_2m + Huma_2m + Ta_10m + Huma_10m + rad_sw_in + rad_sw_out + rad_lw_in + rad_lw_out + Albedo + rad_net + Windspeed_2m + Wind_direction_3m + Windspeed_10m + Wind_direction_10m + PCP,
                data=clim_stat_2)
# Get the imputed values
imputed <-impute.transcan(f, data=clim_stat_2, imputation=1, list.out=TRUE)
# convert the list to the database
imputed.data <- as.data.frame(do.call(cbind,imputed))
imputed.data$date_time_hourly = as.POSIXct(clim_stat_2$date_time_hourly, format = "%Y-%m-%d %H:%M:%S")
# re-check gaps
mice_plot <- aggr(imputed.data, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(imputed.data), cex.axis=.3,
                  gap=3, ylab=c("Missing data","Pattern"))
# 3 - output
#-------------

saveRDS(imputed.data, gw_clim_rdsFN)

# check if sane
clim_stat <- readRDS(gw_clim_rdsFN)
summary(clim_stat)
nrow(clim_stat)


