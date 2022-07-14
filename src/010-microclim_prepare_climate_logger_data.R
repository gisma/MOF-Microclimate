# --------------------------------------------------
# Name: 010_prepare_climate_logger_data
# Author: Lena Perzlmaier, Chris Reudenbach
# Description: load and clean data from data loggers
#              aggregate data into 1h intervals
# Data:   - Ecowitt Daten: (were not included in the model)
#             Beschreibung der Daten: https://nature40.github.io/Nature40DocumentationProject//unit03/04-climate-deployment
#             Datendownload: http://gofile.me/3Z8AJ/OAeFIeLN1
#         - treetalker: 
# Output: "climate_stations_combined.csv" und "climate_stations_combined.rds"
# -------------------------------------------------


# 0 - load packages
#------------------
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# Load the library
appendpackagesToLoad = c("lidR","future","lwgeom","tmap","ggplot2","VIM","Hmisc","mice")
appendProjectDirList =  c("data/lidar/",
                          "data/GRASS/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/",
                          "data/lidar/MOF_lidar_2018/lidar_ctg/",
                          "data/lidar/lidar_orig/",
                          "data/lidar/MOF_lidar_2018/lidar_norm/",
                          "data/lidar/MOF_lidar_2018/lidar_vector/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/pmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/tmetrics/",
                          "data/lidar/MOF_lidar_2018/lidar_raster/topo/",
                          "data/measurements/measure_org",
                          "data/measurements/measure_lev0",
                          "data/measurements/measure_lev1")

source(file.path(root_folder,"src/000-rspatial-setup.R"),echo = TRUE)


gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")
gw_clim_csvFN = file.path(envrmt$path_measure_org, "302_KlimaturmCaldernWiese_complete.csv")

startdate = "2018-01-01"
enddate   = "2022-01-03"

# 3 - start code
#----------------
 

#' ## 3.1 - Ecowitt Dataset
#' #----------------
#' 
#' #Diese Daten bestehen aus 5 Wetterstationen mit je 8 Sensoren (CH01 bis CH08)
#' # 
#' # Es gibt aber eine weitere WeatherStationID: 82EB09C99462A1AEB32109AFC514F264 --> Zu welchem Standort gehört die?
#' 
#' WeatherStationID <- c("E5C37C83C2C1EB2C62BC8FC299D8E281",
#'                       "D3E7B46BC773DA5D0D30F7B94570700F",
#'                       "64B95BE38E36DB75ECE69E193C9CFD62",
#'                       "C90B7356EFDA99C5584B057DA2363EDB",
#'                       "CC354E32668384323E5CB75B8D4A4722",
#'                       "82EB09C99462A1AEB32109AFC514F264")
#' # zugehörige HubID
#' 
#' HubID <- c("GW100A-WIFI243E",
#'            "GW100A-WIFIF774",
#'            "GW100A-WIFIF824",
#'            "GW100A-WIFIFB27",
#'            "GW100A-WIFIFC29",
#'            NA)
#' 
#' # zugehörige Core Study Tree ID
#' cst <- c("cst_00028",
#'          "cst_00002",
#'          "cst_00076",
#'          "cst_00051",
#'          "cst_00050",
#'          NA)
#' 
#' # zugehörige Koordinaten aus der MOFgeoDB (vw_tree_living)
#' x <- c(477706.9,
#'        477754.6,
#'        477770.457,
#'        477830.175,
#'        477818.2,
#'        NA)
#' 
#' y <- c(5632326.1,
#'        5632039.7,
#'        5632254.788,
#'        5632168.856,
#'        5632010.3,
#'        NA)
#' 
#' ID_df <- data.frame("WeatherStationID" = WeatherStationID,
#'                     "HubID" = HubID,
#'                     "cstID" = cst, 
#'                     "X" = x,
#'                     "Y" = y)
#' saveRDS(ID_df, file.path(envrmt$path_measurements, "WeatherStationID.rds"))
#' head(ID_df)
#' 
#' #' 
#' #' Jetzt können die csvs eingelesen werden
#' #' 
#' #' "weather_CC354E32668384323E5CB75B8D4A4722_2020-06-16.csv" [36] und
#' #' "weather_E5C37C83C2C1EB2C62BC8FC299D8E281_2020-09-12.csv" [126] werden nicht eingelesen, da das csv fehlerhaft sind
#' 
#' file_names <- list.files(path = file.path(envrmt$path_measurements, "org/ecowitt"))
#' file_names <- file_names[-c(36,126)]
#' 
#' list <- list()
#' 
#' for (i in 1:length(file_names)) try({
#'   
#'   df_current <- read.table(file = file.path(envrmt$path_measurements, "org/ecowitt/", file_names[i]), 
#'                            header = TRUE, 
#'                            sep = ",", 
#'                            fill = TRUE,
#'                            skipNul = TRUE)
#'   
#'   list[[i]] <- df_current
#' })
#' 
#' 
#' #' Data cleaning:
#' #' 
#' #' * Nur Zeilen die eine gültige PASSKEY haben, werden behalten
#' #' * Nur relevante Spalten werden behalten (siehe colnames)
#' #' * Werte in numerische Werte umwandeln
#' #' * Für Temp muss der Wert größer als 10°F sein (grob geschätzt welche Temperaturen hier im Sommer üblich sind)
#' #' * Alle Zeilen mit NAs werden rausgeschmissen
#' 
#' list2 <- list()
#' 
#' for (i in 1:length(list)){
#'   if (!is.null(list[[i]][["PASSKEY"]]) ){
#'     list2[[i]] <- filter(list[[i]], PASSKEY %in% WeatherStationID)
#'   } 
#' }
#' 
#' colnames <- colnames(list[[49]][1:23])
#' colnames
#' 
#' for (i in 1:length(list2)) {
#'   
#'   list2[[i]] <- list2[[i]][,which(colnames(list2[[i]]) %in% colnames)]
#'   
#'   for (j in 8:23) {
#'     if (!is.null(list2[[i]][[colnames[[j]]]])) {
#'       list2[[i]][[colnames[j]]] <- as.numeric(list2[[i]][[colnames[j]]])
#'     }
#'   }
#'   
#'   for (j in c(8, 10, 12, 14, 16, 18, 20, 22)) {
#'     if (!is.null(list2[[i]][[colnames[[j]]]])) {
#'       list2[[i]] <- filter(list2[[i]], list2[[i]][[colnames[j]]] >= 10)
#'     }
#'   }
#'   
#'   list2[[i]] <- list2[[i]][complete.cases(list2[[i]]),]
#'   
#' }
#' 
#' 
#' #' 
#' #' Dataframes werden zusammengefasst nach PASSKEY
#' 
#' list3 <- vector(mode = "list", length = 6)
#' 
#' for (i in 1:length(list2)) try({
#'   
#'   for (j in 1:length(WeatherStationID)) {
#'     
#'     if (unique(list2[[i]]$PASSKEY) == ID_df$WeatherStationID[j]) {
#'       
#'       list3[[j]] <- plyr::rbind.fill(list3[[j]],list2[[i]])
#'       
#'     }
#'   }
#' })
#' 
#' 
#' 
#' #' 
#' #' Final df
#' 
#' 
#' df <- do.call("rbind.fill", list3)
#' 
#' df_long <- pivot_longer(df, 
#'                         cols = temp1f:humidity8,
#'                         names_to = c("measure","channel"),
#'                         names_pattern   = "(humidity|temp)([1-8])",
#'                         values_to = "measurement",
#'                         values_drop_na = TRUE)
#' 
#' df_ew_final <- pivot_wider(df_long,
#'                            names_from = measure,
#'                            values_from = measurement)
#' 
#' # Umrechnung von Fahrenheit in Celsius
#' df_ew_final$temp <- (df_ew_final$temp -32) * (5/9)
#' 
#' # zu stündlichen Werten aggregieren
#' df_ew_final$date <- as.POSIXct(df_ew_final$dateutc, format = "%Y-%m-%dT%H")
#' 
#' df_ew_final_agg <- df_ew_final %>% dplyr::group_by(date, PASSKEY, channel) %>% 
#'   dplyr::summarise(temp = mean(temp),
#'             humidity = mean(humidity))
#' 
#' # remove passkey "82EB09C99462A1AEB32109AFC514F264" because no tree is associated
#' df_ew_final_agg <- filter(df_ew_final_agg, PASSKEY !="82EB09C99462A1AEB32109AFC514F264")
#' 
#' # Visualize
#' head(df_ew_final_agg)
#' 
#' ggplot(data = df_ew_final_agg, aes(x = date, y = temp, color = PASSKEY, group = channel)) +
#'   geom_point(size = 0.25)
#'  

## 3.2 - Tree Talker Dataset
#----------------

df_tt <- read.table(file.path(envrmt$path_measure_org, "tree_talker/plots.csv"),
                    header = TRUE,
                    sep = ",")
df_tt_f <- df_tt %>% filter(between(year, "2018", "2022"))

sensor_description <- read.table(file.path(envrmt$path_measure_org,"tree_talker/sensor_description.csv"),
                                 header = TRUE,
                                 sep = ",", 
                                 encoding = "UTF-8")

plot_description <- read.table(file.path(envrmt$path_measure_org, "tree_talker/plot_description.csv"),
                               header = TRUE,
                               sep = ",")

cst_tt <- read_sf(file.path(envrmt$path_auxdata, "cst_tt.geojson"))

#' 
#' convert df_tt into long format (with a column defining the plot id and columns with the climatic variables)
#' 
## -------------------------------------------------------------------------------------------------------------------------
df_tt_long <- pivot_longer(df_tt, 
                           cols = tt_air_relative_humidity:ttraw_air_temperature,
                            names_to = c("climate_var"),
                           # names_pattern =  "(tt_[[:alnum:]]*).((?<=\\.).*)",
                            values_to = "measurement",
                            values_drop_na = TRUE)

df_tt_final <- pivot_wider(df_tt_long,
                           names_from = climate_var,
                           values_from = measurement)

head(df_tt_final)

# are the values mean, max and min values always the same?
unique(df_tt_final$tt_air_relative_humidity==df_tt_final$tt_air_relative_humidity_max)
unique(df_tt_final$tt_air_relative_humidity==df_tt_final$tt_air_relative_humidity_min)
unique(df_tt_final$tt_air_temperature ==df_tt_final$tt_air_temperature_max)
unique(df_tt_final$tt_air_temperature ==df_tt_final$tt_air_temperature_min)
# yes! the min/max columns can be deleted

df_tt_final <- df_tt_final[,-c(9, 10, 12, 13)]

df_tt_final$date <- as.POSIXct(df_tt_final$datetime, format = "%Y-%m-%dT%H")

#' it seems that some treetalkers are not correctly calibrated, especially tree_ids: tt_52030105, tt_52030112, tt_52030128, tt_52030129, tt_52030134, tt_52030138
#' Those have values higher than 100 and a break in the data from very high(> 100) to very low (~0).
#' Trees, which have relative humidity values above 100 are set to 100
#' Trees, which have relative humidity values below 15 are set to 15
#' Trees, which have temperaturre values above 55 are set to 55
#' Trees, which have temperaturre values below -25 are set to -25
## -------------------------------------------------------------------------------------------------------------------------


boxplot(df_tt_final$tt_air_temperature)
# there are extreme outlier, which will be rest to boundary values 

df_tt_final = df_tt_final %>%
  mutate(tt_air_temperature = replace(tt_air_temperature,tt_air_temperature>55, 55)) %>%
  as.data.frame()
df_tt_final = df_tt_final %>%
  mutate(tt_air_temperature = replace(tt_air_temperature,tt_air_temperature<-25, 25)) %>%
  as.data.frame()

boxplot(df_tt_final$tt_air_temperature)

ggplot(data = df_tt_final, aes(x = date, y = tt_air_temperature)) +
  geom_point(size = 0.25)+
  facet_wrap(~plotID)
boxplot(df_tt_final$tt_air_relative_humidity)
# there is humidity > 100, which should be reset to 100

df_tt_final = df_tt_final %>%
  mutate(tt_air_relative_humidity = replace(tt_air_relative_humidity,tt_air_relative_humidity>100, 100)) %>%
  as.data.frame()

df_tt_final = df_tt_final %>%
  mutate(tt_air_relative_humidity = replace(tt_air_relative_humidity,tt_air_relative_humidity<15, 15)) %>%
  as.data.frame()

boxplot(df_tt_final$tt_air_relative_humidity)

ggplot(data = df_tt_final, aes(x = date, y = tt_air_relative_humidity)) +
  geom_point(size = 0.25)+
  facet_wrap(~plotID)

ggplot(data = df_tt_final, aes(x = date, y = tt_air_relative_humidity)) +
  geom_point(size = 0.25, na.rm = TRUE)

hist(df_tt_final$tt_air_relative_humidity)
hist(df_tt_final$tt_air_temperature)

# control on gaps
md.pattern(df_tt_final)
mice_plot <- aggr(df_tt_final, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(df_tt_final), cex.axis=.3,
                  gap=3, ylab=c("Missing data","Pattern"))


#' # Combine Datasets
#' 
## -------------------------------------------------------------------------------------------------------------------------
# #str(df_ew_final_agg)
# str(df_tt_final)
# 
# df_final1 <- data.frame(date = df_ew_final_agg$date,
#                         type = "ecowitt",
#                         station_id = df_ew_final_agg$PASSKEY,
#                         channel = df_ew_final_agg$channel,
#                         rel_humidity = df_ew_final_agg$humidity,
#                         temp = df_ew_final_agg$temp)
# 
df_final2 <- data.frame(date = df_tt_final$date,
                        type = "treetalker",
                        station_id = df_tt_final$plotID,
                        rel_humidity = df_tt_final$tt_air_relative_humidity,
                        temp = df_tt_final$tt_air_temperature)
# 
# df_final <- rbind(df_final1, df_final2)

# merge with cst_ids


# trainSites dataset works with tree ID (cst_...)
# need to add tree ID information to trainDF (which currently has the original PASSKEY and serial number as ID)

# ecowitt <- data.frame(station_id = ID_df$WeatherStationID,
#                       cst_id = paste("mof_", ID_df$cstID, sep =""))
# 
# ecowitt[6,2] <- NA

treetalker <- data.frame(station_id = paste("tt_", cst_tt$Serial.Number.Treetalker, sep = ""),
                         cst_id = cst_tt$name)
treetalker <- unique(treetalker)

# tt_52030115 cannot be associated uniquely with on tree ID therefore deleted
treetalker <- filter(treetalker$station_id, station_id != "tt_52030115")

# eco_tt <- rbind(ecowitt, treetalker)

df_final <- merge(df_final2, treetalker, by.x = "station_id", by.y="cst_id", all.x = TRUE, all.y = FALSE)
names(df_final) = c("tree_id","date","type", "rel_humidity", "temp", "ttID")
head(df_final)

md.pattern(df_final)
mice_plot <- aggr(df_final, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(df_tt_final), cex.axis=.3,
                  gap=3, ylab=c("Missing data","Pattern"))


write.csv(df_final, file.path(envrmt$path_measure_lev0, "climate_stations_combined.csv"), row.names = FALSE)
saveRDS(df_final, file.path(envrmt$path_measure_lev0, "climate_stations_combined.rds"))

 
