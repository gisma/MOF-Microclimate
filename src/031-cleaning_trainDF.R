#
# Name: 040_clean_dataset
# Author: Lena Perzlmaier, Chris Reudenbach
# Description:  creates the training data frame
#             - stack topo and forest structure  -> pred_stack
#             - get coordinates of the trainSites and extract information of predStack
#             - merge information of the step before with data from the climate station
#             - merge with radiation information
# Input:   see source files
# 
#


# ---- project setup ----
# ----project setup----
require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("ggplot2","VIM","Hmisc","mice","lubridate")

# define additional subfolders comment if not needed

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

# MANDANTORY: calling the setup script also DO NOT change this line
source(file.path(envimaR::alternativeEnvi(root_folder = root_folder),"src/000-rspatial-setup.R"))

gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")
gw_clim_csvFN = file.path(envrmt$path_measure_org, "302_KlimaturmCaldernWiese_complete.csv")
combo_clim_rdsFN = file.path(envrmt$path_measure_lev0, "climate_stations_combined.rds")

# load trainDF from script 30
trainDF <- readRDS(file.path(envrmt$path_measure_lev0, "trainDFmc.rds"))


# # have a look at ecowitt data
# # it is only available for short time periods
c
# only keep treetalker
trainDFclean <- filter(trainDF, c(type == "treetalker" | tree_id == "gw_station" ))
trainDFclean$type = NULL
x = trainDFclean$x
y = trainDFclean$y

##---  cleaning the training data for a random forest model training
traintmp = trainDFclean[ , !(names(trainDFclean) %in% c("cst_id.y","zq5","WSDiag_2","WSDiag_1","cell","doy_hour_id","Albedo","date","cst_id.x","station_id","channel","n","area","x","y","cell_coverage_fraction","heatflux_soil","cst_id.y"))]
summary(traintmp)

## ---- 3.x cloudiness index ----

# calculate the ratio of observed versus calculate irradiance to obtain an indicator of cloudiness

cloud <- vector(length = nrow(traintmp))
for (i in 1:nrow(traintmp)) {
  if(traintmp[i,]$rad <= 0){
    cloud[i] <- 0
  } else {
    cloud[i] <- traintmp[i,]$rad_sw_in/traintmp[i,]$rad
  }
}
summary(cloud)
traintmp$ttID=NULL
hist(cloud,breaks="Scott")
traintmp$cloudiness <- cloud
traintmp <- traintmp %>% mutate(day = if_else(rad > 0, 1, 0))
cness = traintmp$cloudiness
cloudi_daily <- traintmp %>%  group_by(doy) %>% summarise(cloudi_daily = mean(cloudiness, na.rm = TRUE))
plot(cloudi_daily)




##--- filter zero or near-zero values
# nzv = nearZeroVar(traintmp)
# names(traintmp[,nzv])
# if (length(nzv) > 0) traintmp = traintmp[, -nzv]

##--- filter correlations that are > cor_cutoff
#filt = findCorrelation(cor(subset(traintmp, select=-c(tree_id)), use = "complete"), cutoff = 0.98)
#traintmp = traintmp[,-filt]

##--- re-add the necessary variables for model training

trainDFclean$x = x
trainDFclean$y = y


##--- remove rows with NA
traintmp = traintmp[complete.cases(traintmp) ,]

##--- check manually if there are still NA values around
summary(traintmp)
sapply(traintmp, function(y) sum(length(which(is.na(y)))))

##--- define  the classes to be factor - is obligate for rr
traintmp$tree_id = as.factor(traintmp$tree_id)

# the variables rad_Klimastation and cloudiness have "spikes" in their distribution
hist(traintmp$rad, 24)
hist(traintmp$cloudiness, 24)

# linear mixed models have issues with different scales
# lmer function gives warning, that some predictors have very different scales
# rf is quite robust to different scales, but svm need scaling

preProcValues <- caret::preProcess(traintmp, method = c("center", "scale"))
traintmpscaled = predict(preProcValues, traintmp)

## -----4 output------

saveRDS(traintmp, file.path(envrmt$path_measure_lev1, "trainDFmcClean.rds"))
saveRDS(traintmpscaled, file.path(envrmt$path_measure_lev1, "trainDFmcScaled.rds"))
saveRDS(preProcValues, file.path(envrmt$path_measure_lev1, "preProcValues.rds"))


## ----5 visualizations -----

trainDF <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcClean.rds"))

ggplot(trainDF, aes(x = rad, y = temp))+
  geom_point( aes(color = hour))+
  facet_wrap(vars(tree_id))

insane_tt = c( "mof_cst_00002", "mof_cst_00003" , "mof_cst_00010" , "mof_cst_00019",
             "mof_cst_00021" , "mof_cst_00027" , "mof_cst_00030" , "mof_cst_00033",
             "mof_cst_00037", "mof_cst_00042", "mof_cst_00043",  "mof_cst_00009",
             "mof_cst_00025", "mof_cst_00023", "mof_cst_00044",  "mof_cst_00018",
             "mof_cst_00050", "mof_cst_00052", "mof_cst_00076", "mof_cst_00049",
             "mof_cst_00040", "mof_cst_00045") 


t=trainDF[!trainDF$tree_id %in% insane_tt, ]
ggplot(t, aes(x = rad, y = temp))+
  geom_point( aes(color = hour))+
  facet_wrap(vars(tree_id))
saveRDS(t, file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))
  