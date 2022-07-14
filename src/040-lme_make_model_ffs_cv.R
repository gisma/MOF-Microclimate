#------------------------------------------------------------------------------
# Name: 39_forward_feature_selection
# Author: Lena Perzlmaier
# Description:  find the best model, with forward feature selection using AIC, BIC and 
#               k fold crossvalidation (with RMSE) as performance measures
# 
#------------------------------------------------------------------------------

# 0 - load packages
#------------------

require(envimaR)
require(rprojroot)
root_folder = find_rstudio_root_file()

# define  additional packages comment if not needed
appendpackagesToLoad = c("lme4", "nlme", "CAST", "ggplot2", "lattice")

# define additional subfolders comment if not needed
tmpPath = paste0("data/lidar/",prefix,"/") 
appendProjectDirList =  c("data/lidar/",
                          "data/GRASS/",
                          "data/data_lev1/models/",
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

predstack_forest_metricsFN = file.path(envrmt$path_lidar_raster, "pred_forest_structure.tif")
gw_clim_rdsFN = file.path(envrmt$path_measure_lev1, "klimastation_wiese_hourly.rds")

forest_structure_metrics <- rast(predstack_forest_metricsFN)  
clim_stat <- readRDS(gw_clim_rdsFN)
trainDFclean <- readRDS(file.path(envrmt$path_measure_lev1, "trainDFmcCleanSANE.rds"))
names(trainDFclean)

traintmp = trainDFclean
##--- filter correlations that are > cor_cutoff to reduce higly correlated params and cpu time
filt = findCorrelation(cor(subset(traintmp, select=-c(tree_id,cst_id,doy_hour_id)), use = "complete"), cutoff = 0.9)
names(traintmp[,filt])
traintmp = traintmp[,-filt]
#traintmp$Ta_10m = trainDFclean$Ta_10m
#traintmp$Ta_2m = trainDFclean$Ta_2m
traintmp$LAD_metrics_MOF = trainDFclean$LAD_metrics_MOF
traintmp$tree_cluster = trainDFclean$tree_cluster
traintmp$entropy_metrics_MOF = trainDFclean$entropy_metrics_MOF
traintmp$dtm = trainDFclean$dtm
predictors_all = c(names(traintmp[!names(traintmp) %in% c("rel_humidity","tree_id","hour","doy","cst_id","doy_hour_id")]),"chm * doy", "chm * hour", "tree_cluster * doy","tree_cluster * hour" )

# 3 - start code 
# ----------------

# --- 3.1 predictor sets ----

predictors_interaction <- c("Huma_10m", "Huma_2m", "Ta_10m", "aspect","aspect * doy", "aspect * hour",  "slope", "slope * doy", "slope * hour","rad","rad * hour","rad * doy",
                            "chm", "chm * doy", "chm * hour", "tree_cluster", "tree_cluster * doy","tree_cluster * hour" ,"zpcum2","zpcum2 * doy","zpcum2 * hour" )



# --- forward feature selection ----

fwd_fs_RMSE_interaction_hour_huma <- ffs_cv_rmse(traintmp,
                                           dep = "rel_humidity",
                                           sp_var= "tree_id",
                                           t_var = "hour",
                                           vars = predictors_interaction,
                                           random = "(1|doy/hour)",
                                           nfolds = 5)

#plot(traintmp[traintmp$doy == 125,]$Ta_2m,traintmp[traintmp$doy == 125,]$rel_humidity)

fwd_fs_RMSE_interaction_doy_huma <- ffs_cv_rmse(traintmp,
                                            dep = "rel_humidity",
                                            sp_var= "tree_id",
                                            t_var = "doy",
                                            rmse_th = 0.05,
                                            vars = predictors_interaction,
                                            random = "(1|doy/hour)",
                                            nfolds = 5)
saveRDS(fwd_fs_RMSE_interaction_hour_huma,paste0(envrmt$path_models,"fwd_fs_RMSE_interaction_hour_huma.rds"))
saveRDS(fwd_fs_RMSE_interaction_doy_huma,paste0(envrmt$path_models,"fwd_fs_RMSE_interaction_doy_huma.rds"))


fwd_fs_RMSE_all_hour_huma <- ffs_cv_rmse(traintmp,
                                   dep = "rel_humidity",
                                   sp_var= "tree_id",
                                   t_var = "hour",
                                   rmse_th = 0.05,
                                   vars = predictors_all,
                                   random = "(1|doy/hour)",
                                   nfolds = 5)
saveRDS(fwd_fs_RMSE_all_hour_huma,paste0(envrmt$path_models,"fwd_fs_RMSE_all_hour_huma.rds"))



# performance c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# c("Ta_10m", "rad", "rad_lw_out", "chm * hour", "Windspeed_10m", "rad_lw_in", "Ta_2m", "p2th", "rel_humidity", "rad_sw_out")
# c(2.71647443936501, 2.43646042480053, 2.2438790223352, 2.2035154715153, 2.17804019979317, 2.16392323967337, 2.15370572630087, 2.1518178330618, 2.15006951321657, 2.14865751987541)

# fwd_fs_RMSE_all_doy <- ffs_cv_rmse(traintmp,
#                                    dep = "temp",
#                                    sp_var= "tree_id",
#                                    t_var = "doy",
#                                    vars = predictors_all,
#                                    random = "(1|doy/hour)",
#                                    nfolds = 5)
# 

# saveRDS(fwd_fs_RMSE_all_hour,paste0(envrmt$path_models,"fwd_fs_RMSE_all_hour.rds"))


# --- 3.3 resulting models ----

# formula <- paste("temp", " ~ ", 
#                  paste(predictors_interaction, collapse = " + "), 
#                        "+ (1|doy/hour)")
# lme <- lmer(formula, data = train)

interaction_lme <- lmer(temp ~ Ta_10m + rad*hour +  rad*doy + zpcum2*hour + tree_cluster*hour + aspect +(1|doy/hour), traintmp)
summary(interaction_lme)
saveRDS(interaction_lme, file.path(envrmt$path_models, "/interaction_lme.rds"))
BIC(interaction_lme)


all_lme = lmer(temp  ~  Ta_10m + rad + rad_lw_out + chm + hour + Windspeed_10m + rad_lw_in + p2th + (1|doy/hour), trainDFclean)
summary(all_lme)
saveRDS(all_lme, file.path(envrmt$path_models, "/all_lme.rds"))
summary(all_lme)
BIC(all_lme)


