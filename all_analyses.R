rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(raster)
library(ggplot2)
library(dplyr)
library(tidymodels)
library(tidyverse)
library(dismo)
library(rgbif)
library(tmap)
library(CoordinateCleaner)
library(stringi)
library(phytools)
library(rnaturalearth)
library(patchwork)
library(viridis)
library(correlation)
library(phylosignal)
library(phylobase)
library(caper)

data <- read.csv("data/MasterTable11.csv", strip.white = TRUE, 
                 na.strings = c("", "NA"), stringsAsFactors = TRUE)
data$KBTribe[data$KBGenus_Sp == "Rhodnius ecuadoriensis"] <- "Rhodniini"

# function to fix names
fix_species_name <- function(column) {
  require(stringr)
  column_new <- stringr::str_replace(column, "(?s) .*", "sp")
  column_new
}

# host
data$Host_Genus_sp <- str_replace(as.character(data$Host_Genus_sp), "NA NA", 'NA')
data$Host_Genus_sp <- str_replace(as.character(data$Host_Genus_sp), "(?s) NA", " sp")
data$Host_Genus_sp[which(data$Host_Genus_sp == "NA rattus")] <- "Rattus rattus"
data$Host_Genus_sp[which(data$Host_Genus_sp == "Felis spp.")] <- "Felis sp"
data$Host_Genus_sp <- as.factor(data$Host_Genus_sp)
data$Host_Genus_sp[which(data$Host_Genus_sp == "NA")] <- NA
data$Host_genus[which(is.na(data$Host_genus) & !is.na(data$Host_Genus_sp))] <- "Rattus"

# bugs
data$KBGenus_Sp <- str_replace(as.character(data$KBGenus_Sp), "NA NA", 'NA')
data$KBGenus_Sp <- str_replace(as.character(data$KBGenus_Sp), "(?s) NA", " sp")
data$KBGenus_Sp <- as.factor(data$KBGenus_Sp)

# country 
data$Country <- str_replace(as.character(data$Country), "Panamá", "Panama")
data$Country <- str_replace(as.character(data$Country), "Guyana Francesa", "French Guiana")

# estimating species-level and genus-level DBR for each species
# with NAs 
dietbreadthdt <- data %>% 
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR = n_distinct(Family, na.rm = FALSE),
            orderDBR = n_distinct(Order, na.rm = FALSE),
            classDBR = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup()

# without NAs
dietbreadthdtnoNA <- data %>% 
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR = n_distinct(Host_Genus_sp, na.rm = TRUE),
            genusDBR = n_distinct(Host_genus, na.rm = TRUE),
            familyDBR = n_distinct(Family, na.rm = TRUE),
            orderDBR = n_distinct(Order, na.rm = TRUE),
            classDBR = n_distinct(Class, na.rm = TRUE)) %>%
  ungroup()

# cleaning to remove T. rubrofasciata that is wrong 
dietbreadthdt_clean <- dietbreadthdt %>% ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR,
                genusDBR,
                familyDBR,
                orderDBR,
                classDBR) %>%
  na.omit() #%>% filter(species != "Triatoma rubrofasciata")

dietbreadthdtnoNA_clean <- dietbreadthdtnoNA %>% ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR,
                genusDBR,
                familyDBR,
                orderDBR,
                classDBR) %>%
  na.omit() #%>% filter(species != "Triatoma rubrofasciata")

genus_name_unique <- as.list(unique(sub('(^\\w+)\\s.+', '\\1', dietbreadthdt$KBGenus_Sp)))
species_name_unique <- as.list(unique(sub('(^\\w+)\\s(+)', '\\2', dietbreadthdt$KBGenus_Sp)))
species_name_unique <- unique(levels(dietbreadthdt$KBGenus_Sp))

gbif_taxon_keys <- matrix(ncol = 2, nrow = length(species_name_unique))
gbif_taxon_keys[, 1] <- species_name_unique
for (i in 1:nrow(gbif_taxon_keys)) {
  try(
   gbif_taxon_keys[i, 2] <- name_backbone(gbif_taxon_keys[i, 1])$speciesKey
  )
}
gbif_taxon_keys <- gbif_taxon_keys[complete.cases(gbif_taxon_keys[, 2]), ]

gbif_downloads <- c()
syn <- list()
for (i in 1:nrow(gbif_taxon_keys)) {
  gbif_data <- occ_download(pred("taxonKey", gbif_taxon_keys[i, 2]),
                                 pred("hasGeospatialIssue", FALSE),
                                 pred("hasCoordinate", TRUE),
                                 pred_not(pred_in("basisOfRecord", "FOSSIL_SPECIMEN")),
                                 format = "SIMPLE_CSV")

  occ_download_wait(gbif_data[1])

  gbif_download <- gbif_data %>%
    occ_download_get() %>%
    occ_download_import() %>%
    setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
    cc_sea() %>% # remove from ocean 
    cc_val() %>%
    cc_zero() %>%
    glimpse()

  gbif_data_clean <- gbif_download %>% 
    dplyr::select(longitude = decimallongitude, latitude = decimallatitude, 
                  genus, species, country = countrycode, occurrencestatus,
                  basisofrecord, coordinateprecision, 
                  coordinateuncertaintyinmeters)

  syn[[i]] <- unique(gbif_data_clean$species)

  gbif_data_clean$species <- gbif_taxon_keys[i, 1]

  gbif_downloads <- bind_rows(gbif_downloads, gbif_data_clean)

}
#save(gbif_downloads, file = "data/gbif_downloads.RData")

# subsetting to only include species which we have the diet breadth
gbif_data_clean2 <- gbif_downloads[gbif_downloads$species %in% species_name_unique, ]

# subsetting to only include species in the countries of America's
gbif_data_clean2 <- gbif_data_clean2[!(gbif_data_clean2$country %in% c("CN", "ID", "IN", "LK", "PH", "VN")), ]

# selecting the variables 
gbif_data_clean2 <- gbif_data_clean2 %>% dplyr::select(longitude, latitude, 
  genus, species, country)

# need to remove all NA's or function will not work
gbif_data_clean2 <- na.omit(gbif_data_clean2) %>% droplevels(.)

write.csv(gbif_data_clean2, "data/Full_GBIF_data.csv")
#gbif_data_clean2 <- read.csv("data/Full_GBIF_data.csv") 

# species into polygons
# query current environmental data for our location
clim.current <- raster::getData("worldclim", var = "bio", res = 10)

# raster with the same dimensions as above
r <- raster(nrows = 900, ncols = 2160, ymn = -60, vals = 1:(2160*900))

# SDM modelling 
sdm_rf_manual <- function(data, 
                          mtry = c(3, 4, 5), 
                          trees = c(100, 300, 500), 
                          threshold = 0.2,
                          split_prop = 0.75) {
  outputlist <- list()
  # defining the extent of latin america + USA
  LatAm.ext <- raster::extent(-125, -28, -70, 48) 
  geo_data2 <- data
  # selecting coordinates 
  sp.occ <- cbind.data.frame(geo_data2$longitude, geo_data2$latitude)
  
  # creating presence/absence dataframe 
  presence <- cbind(sp.occ, status = 1)
  colnames(presence)[1:3] <- c("lon", "lat", "status")
  
  absence <- cbind(randomPoints(LatAm_clim, (nrow(presence)*10)), status = 0)
  colnames(absence)[1:3] <- c("lon", "lat", "status")
  
  presabs <- rbind(presence, absence)
  
  # extracting environmental variables in df format
  climdt <- cbind(presabs[3], raster::extract(LatAm_clim, presabs[-3]))
  climdt$status <- as.factor(climdt$status)
  climdt <- na.omit(climdt)
  
  # fitting the ML model
  # split the data into trainng (75%) and testing (25%)
  library(tidymodels)
  climdt_split  <- initial_split(climdt, prop = split_prop)
  climdt_split
  
  # extract training and testing sets
  climdt_train <- training(climdt_split)
  climdt_test <- testing(climdt_split)
  
  # create CV object from training data
  climdt_cv <- vfold_cv(climdt_train)
  
  # define the recipe, which consists of the formula (outcome ~ predictors)
  rec <- recipe(status ~ bio1 + bio4 + bio7 + bio12 + bio15,
                data = climdt_train)

  climdt_recipe <- rec %>%
    step_impute_knn(all_predictors()) %>%
    # and some pre-processing steps
    step_normalize(all_numeric())
  
  climdt_train_preprocessed <- climdt_recipe %>%
    # apply the recipe to the training data
    prep(climdt_train) %>%
    # extract the pre-processed training dataset
    juice()
  
  # specify the model 
  rf_model <- 
    # specify that the model is a random forest
    rand_forest() %>%
    # specify that the `mtry` parameter needs to be tuned
    set_args(mtry = tune(), trees = tune()) %>%
    # select the engine/package that underlies the model
    set_engine("ranger", importance = "impurity") %>%
    # choose either the continuous regression or binary classification mode
    set_mode("classification") 
  
  # set the workflow
  rf_workflow <- workflow() %>%
    # add the recipe
    add_recipe(climdt_recipe) %>%
    # add the model
    add_model(rf_model)
  
  # specify which values want to try
  rf_grid <- expand.grid(mtry = mtry, trees = trees)
  # extract results
  rf_tune_results <- rf_workflow %>%
    tune_grid(resamples = climdt_cv, # CV object
              grid = rf_grid, # grid of values to try
              metrics = metric_set(accuracy, roc_auc) # metrics we care about
              )
  # print results
  rf_tune_results %>%
    collect_metrics()
  
  # extract the best values 
  param_final <- rf_tune_results %>%
    select_best(metric = "roc_auc")
  
  # finalise workflow
  rf_workflow <- rf_workflow %>%
    finalize_workflow(param_final)
  
  # evaluating the model
  rf_fit <- rf_workflow %>%
    # fit on the training set and evaluate on test set
    last_fit(climdt_split)
  
  test_performance <- rf_fit %>% collect_metrics()
  
  # generate predictions from the test set
  test_predictions <- rf_fit %>% collect_predictions()
  
  # fitting and using the final model 
  final_model <- fit(rf_workflow, climdt) # full dataset used for the split
  
  # predicting across the entire Latin America + USA
  # Latin America climate data (entire region for prediction)
  sdmpred_clim <- na.omit(data.frame(raster::extract(LatAm_clim, LatAm.ext)))
  coord_clim <- data.frame(coordinates(LatAm_clim),
                           ifelse(!is.na(raster::extract(LatAm_clim, LatAm.ext)) == FALSE, NA, 1)) %>%
    na.omit() %>% select(x, y)
  
  # prediction
  species_dist <- data.frame(coord_clim, 
                             dist = predict(final_model, 
                                            new_data = sdmpred_clim,
                                            type = "prob"))  
  colnames(species_dist) <- (c("x", "y", ".pred_0", ".pred_1"))
  
  species_dist_final <- species_dist %>% filter(.pred_1 >= threshold)
  
  # creating Spatial objects to return
  spatial_presence <- SpatialPoints(coords = cbind(presence$lon, presence$lat),
                                    proj4string = CRS("+init=epsg:4326"))
  
  # creating a spatial training set of presence
  sdm_spatial <- SpatialPoints(coords = cbind(species_dist_final$x, 
                                              species_dist_final$y),
                               proj4string = CRS("+init=epsg:4326"))
  sdm_spatial$pred_prob <- species_dist_final$.pred_1
  
  outputlist <- list(input = sp.occ,
                     final_model = final_model,
                     model_performance = test_performance,
                     spatial_output_presence = sdm_spatial,
                     spatial_observed_data = spatial_presence)
  
  return(outputlist)

}

# DO NOT RUN! too long (1.438433 hours) load data above:
sdm_model_fulldt <- gbif_data_clean2 %>% 
  select(-country) %>% 
  nest(data = -species) %>%
  mutate(sdm = purrr::map(data, ~sdm_rf_manual(.x))) 
    #There were issues with some computations   A: x7
    #There were issues with some computations   A: x5
    #There were issues with some computations   A: x3
    #There were issues with some computations   A: x1
    #There were issues with some computations   A: x1

#save(sdm_model_fulldt, file = "data/sdm_model_fulldt.RData")
#load(file = "data/sdm_model_fulldt.RData")

calculate_SDM_area <- function(data, threshold = 0.5) {
  unprojected_obj <- rasterFromXYZ(data)
  crs(unprojected_obj) <- projection(data)
  area_proj <- raster::area(unprojected_obj)
  avg_area_proj <- cellStats(unprojected_obj > threshold, sum) * area_proj
  output <- list(min_value = area_proj@data@min,
                 max_value = area_proj@data@max,
                 average_area_km2 = mean(getValues(avg_area_proj)))
  return(output)
}

# estimating niche distribution based on SDM
area_model_fulldt <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>%  
  unnest(sdm) %>%
  group_by(species) %>% 
  slice(4) %>%
  mutate(area_sdm = purrr::map(sdm, ~calculate_SDM_area(.x))) %>% 
  dplyr::select(-sdm) %>%
  unnest(area_sdm) %>% 
  slice(3) %>% 
  unnest(area_sdm)

# getting the parameters from the model fit 
param_random_forest_fulldt <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>% 
  unnest(sdm) %>% 
  group_by(species) %>% 
  slice(3) %>% 
  unnest() %>% 
  dplyr::select(-.estimator) %>% 
  tidyr::spread(., '.metric', '.estimate')

# plotting to check sdm

gbif_data_clean2_sf <- sf::st_as_sf(gbif_data_clean2,
                                    coords = c("longitude", "latitude"),
                                    crs = "+proj=longlat +datum=WGS84")

raster_obj <- function(data) {
  raster_unnest <- rasterFromXYZ(data)
  crs(raster_unnest) <- projection(data)
  return(raster_unnest)
}

raster_obj <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>% 
  unnest(sdm) %>% 
  group_by(species) %>% 
  slice(4) %>% 
  ungroup() %>% 
  mutate(raster_obj = purrr::map(sdm, ~raster_obj(.x))) %>% 
  dplyr::select(-sdm)

points_obj <- raster_obj %>%
  mutate(points_obj = purrr::map(raster_obj, ~rasterToPoints(.x))) %>% 
  dplyr::select(-raster_obj) %>%
  mutate(points_obj_df = purrr::map(points_obj, ~data.frame(.x))) %>% 
  dplyr::select(-points_obj)


worldmap <- ne_countries(scale = "medium", type = "map_units",
                         returnclass = "sf", continent = c("south america", 
                                                           "north america",
                                                           "central america"))

pdf("figures/ras.pdf", width = 12, height = 6)
for (i in 1:length(unique(gbif_data_clean2_sf$species))) {
  spp <- unique(gbif_data_clean2_sf$species)[i]

  occ <- ggplot() +
    geom_sf(data = worldmap, color = "gray", fill = "lightgray") +
    xlim(-180, 0) + 
    geom_sf(data = gbif_data_clean2_sf %>% filter (species == spp)) + 
    theme_void() +
    labs(title = spp)

  ras <- ggplot() +
    geom_sf(data = worldmap, color = "gray", fill = "lightgray") +
    geom_raster(data = ((points_obj %>% filter (species == spp) %>% select(points_obj_df))$points_obj_df[[1]]), 
                aes(x = x, y = y, fill = layer)) +
    xlim(-180, 0) + 
    scale_fill_viridis(name = "pred_prob") + 
    theme_void() 

  final <- wrap_plots(occ, ras)
  print(final)
}
dev.off()

# diet breadth data 

data <- data %>%
  filter(KBGenus_Sp %in% sdm_model_fulldt$species) %>%
  droplevels()

dietbreadthdt <- dietbreadthdt %>%
  filter(KBGenus_Sp %in% sdm_model_fulldt$species) %>%
  droplevels()

dietbreadthdtnoNA <- dietbreadthdtnoNA %>%
  filter(KBGenus_Sp %in% sdm_model_fulldt$species) %>%
  droplevels()

dietbreadthdt_clean <- dietbreadthdt_clean %>%
  filter(species %in% sdm_model_fulldt$species) %>%
  droplevels()

dietbreadthdtnoNA_clean <- dietbreadthdtnoNA_clean %>%
  filter(species %in% sdm_model_fulldt$species) %>%
  droplevels()

# join area with DBR
areaDBR_fulldt <- area_model_fulldt %>% inner_join(., dietbreadthdt_clean, by = 'species') 

# raster with the same dimensions as above
r <- raster(nrows = 900, ncols = 2160, ymn = -60, vals = 1:(2160*900))

# annual actual evapotranspiration
aet <- raster("data/aet2.tif") # 21600 x 43200
aet_ras <- resample(aet, r)
#save(aet_ras, file = "data/aet_ras_resampled.RData")
#load(file = "data/aet_ras_resampled.RData")

# annual aridity index
ai <- raster("data/ai2.tif")
ai_ras <- resample(ai, r)
#save(ai_ras, file = "data/ai_ras_resampled.RData")
#load(file = "data/ai_ras_resampled.RData")

# Priestley-Taylor alpha coefficient
alpha <- raster("data/alpha2.tif")
alpha_ras <- resample(alpha, r)
#save(alpha_ras, file = "data/alpha_ras_resampled.RData")
#load(file = "data/alpha_ras_resampled.RData")

# annual potential evapotranspiration
etyr <- raster("data/etyr2.tif")
etyr_ras <- resample(etyr, r)
#save(etyr_ras, file = "data/etyr_ras_resampled.RData")
#load(file = "data/etyr_ras_resampled.RData")

# net primary productivity 
npp <- raster("data/npp08.18.tif")
npp_ras <- resample(npp, r)
#save(npp_ras, file = "data/npp_ras_resampled.RData")
#load(file = "data/npp_ras_resampled.RData")

clim.current2 <- stack(clim.current, aet_ras, ai_ras, alpha_ras, etyr_ras,
                       npp_ras)
env_vals <- raster::extract(clim.current2, 1:length(r[]))  

# predicted niche 

# calculate both realised and predicted niche breadth 
  # (former: spatial points, latter: SDM predictions)
estimate_nichebreadth <- function(.x, threshold = 0.5, sp_clim = env_vals) {

  if ((class(.x) == "SpatialPointsDataFrame") == TRUE) {

    data <- rasterFromXYZ(.x)
    obj_crs <- crs(data)

    # setting up cells as NA below threshold
    data[data < threshold] <- NA
  
  } else if((class(.x) == "SpatialPoints") == TRUE) {
  
    data <- .x
  
  } else {
    warning("Input needs to be either a SpatialPoints or SpatialPointsDataFrame")
  }

  # extracting climate data
  r <- raster(ncols = 2160, nrows = 900, ymn = -60, vals = rep(1, 2160*900))

  raster_i <- resample(data, r)

  rastercells <- which(getValues(raster_i) > 0)
  
  values_env <- sp_clim[rastercells, ]
  values_env <- as.data.frame(values_env)

  env_avg <- numeric()
  for (i in 1:ncol(values_env)) {

    if (nrow(values_env) == 1) {
      ifelse(is.na(values_env[i]),
             env_avg[i] <- NA,
             env_avg[i] <- values_env[i])
    } else {
      ifelse(all(is.na(values_env[, i])),
             env_avg[i] <- NA,
             env_avg[i] <- mean(values_env[, i], na.rm = TRUE))
    }

  }
  names(env_avg) <- colnames(values_env)

  temp_niche_breadth <- max(values_env$bio5, na.rm = T) - 
    min(values_env$bio6, na.rm = T)
  names(temp_niche_breadth) <- "temp_niche_breadth"
  prec_niche_breadth <- max(values_env$bio16, na.rm = T) - 
    min(values_env$bio17, na.rm = T)
  names(prec_niche_breadth) <- "prec_niche_breadth"

  return(c(env_avg, temp_niche_breadth, prec_niche_breadth))
}

predicted_nichebreadth_fulldt <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>%  
  unnest(sdm) %>%
  group_by(species) %>% 
  slice(4) %>% 
  mutate(pred_nichebreadth = purrr::map(sdm, ~estimate_nichebreadth(.x))) %>%
  dplyr::select(-sdm) %>%
  unnest_wider(pred_nichebreadth)

# putting all info into one dataset
area_niche_DBR_fulldt <- areaDBR_fulldt %>% 
  inner_join(., predicted_nichebreadth_fulldt, by = "species")

write.csv(area_niche_DBR_fulldt, "data/final_data_25oct23.csv", row.names = F)

# laoding again in a different format
area_niche_DBR_fulldt <- read.csv("data/final_data_25oct23.csv")

# loading phylogenies (MCC and 1000 random trees from posterior distribution)
tree <- read.nexus("data/tr_cal_f.tre")
trees <- read.nexus("data/1k_random_tria.trees")

# eliminate the samples not in the phylogeny and vice-versa
area_niche_DBR_fulldt$species <- 
  stri_replace_all_fixed(area_niche_DBR_fulldt$species, " ", "_")

area_niche_DBR_fulldt_clean <- area_niche_DBR_fulldt[area_niche_DBR_fulldt$species %in% tree$tip.label, ]

tree <- drop.tip(tree, setdiff(tree$tip.label, area_niche_DBR_fulldt_clean$species))
trees.pruned <- list()
for (i in 1:1000) {
  trees.pruned[[i]] <- drop.tip(trees[[i]], setdiff(trees[[i]]$tip.label,
                                                    area_niche_DBR_fulldt_clean$species))
}

# calculating percentages of each host
dietbreadth_perc <- data %>% 
  filter(KBGenus_Sp %in% stri_replace_all_fixed(area_niche_DBR_fulldt_clean$species, "_", " ")) %>%
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(species = KBGenus_Sp, host_family = Family) %>% 
  summarise(feeds = sum(N_feeds, na.rm = TRUE)) %>%
  filter(species %in% area_model_fulldt$species) %>% 
  droplevels()

pdf("figures/new/Fig1.pdf", height = 10)
dietbreadth_perc %>% 
  filter(!is.na(host_family)) %>%
  ggplot(aes(x = species, y = feeds, fill = host_family)) + ylab("Proportion of diet") +
  geom_bar(position = "fill", stat = "identity") + 
  coord_flip() + 
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  guides(fill = FALSE) + 
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "none") + 
  #theme(legend.position = "right", legend.direction = "vertical")
  geom_text(aes(y = feeds, label = toupper(as.character(host_family))), size = 2, position = position_fill(vjust = 0.5), color = "white")
dev.off()

######################## FIGURE PHYLOGENY

pdf("figures/Fig1_scale.pdf")
revts(ggtree(tree, layout = "fan", open.angle = 180, ladderize = F) + geom_tiplab(size = 1.6)) +
  coord_geo_polar(dat = "epochs", lwd = NULL) +
  scale_x_continuous(breaks = seq(-60, 0, 10), labels = abs(seq(-60, 0, 10)))
dev.off()

DBR_area <- dat[, c("speciesDBR", "area")]
DBR_area <- DBR_area[tree$tip.label, ]
DBR_area_norm <- DBR_area
DBR_area_norm$speciesDBR <- DBR_area_norm$speciesDBR/max(DBR_area_norm$speciesDBR)
DBR_area_norm$area <- DBR_area_norm$area/max(DBR_area_norm$area)

h <- max(nodeHeights(tree))
m <- ncol(DBR_area)
cols <- wes_palette("FantasticFox1")[c(2, 4)]
xlim <- ylim <- 1.4 * c(-h, h) + c(-1, 1) * 0.1 * m * h + 0.02 * c(-h, h)
ylim <- c(0, ylim[2])

pdf("figures/Fig1_tree.pdf", width = 10)
plotTree(tree, type = "fan", lwd = 1, xlim = xlim, ylim = ylim, ftype = "i",
         fsize = 0.5, part = 0.5)
for (i in 1:ncol(DBR_area_norm)) {
    tt <- tree
    tt$edge.length[which(tt$edge[, 2] <= Ntip(tt))] <-
      tt$edge.length[which(tt$edge[, 2] <= Ntip(tt))] +
        0.47 * h + (i - 1) * 0.1 * h
    plotTree(tt, color = "transparent", type = "fan", xlim = xlim, ylim = ylim,
             lwd = 1, ftype = "off", add = TRUE, part = 0.5)
    
    pp1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    tt$edge.length[which(tt$edge[, 2] <= Ntip(tt))] <-
      tt$edge.length[which(tt$edge[, 2] <= Ntip(tt))] + 0.09 * h

    #plotrix::draw.circle(0, 0, radius = h + 0.25 * h + (i - 1) * 0.1 * h,
                         #border="#505050")
    plotTree(tt, color = "transparent", type = "fan", xlim = xlim, ylim = ylim,
             lwd = 1, ftype = "off", add = TRUE, part = 0.5)
    
    pp2 <-get("last_plot.phylo", envir = .PlotPhyloEnv)
    par(lend = 1)
    for(j in 1:Ntip(tree)){
        ii <- which(rownames(DBR_area_norm) == tree$tip.label[j])
        dx <- (pp2$xx[j] - pp1$xx[j]) * DBR_area_norm[ii, i]
        dy <- (pp2$yy[j] - pp1$yy[j]) * DBR_area_norm[ii, i]
        #lines(pp1$xx[j]+c(0,1.05*dx),pp1$yy[j]+c(0,1.05*dy),lwd=10,
        #   col=par()$fg)
        lines(pp1$xx[j] + c(0, dx), pp1$yy[j] + c(0, dy), lwd = 10, 
              col = cols[i])
    }
}

xx <- rep(0.65 * par()$usr[4], m)
yy <- 0.95 * par()$usr[4] - 1.5 * 0:(m-1) * strheight("W")
text(xx, yy, c("DBR", "Range area"), pos = 4, cex = 0.8)
scale.bar <- apply(DBR_area_norm, 2, max, na.rm = TRUE) * 0.09 * h
xx2 <- xx - scale.bar
segments(x0 = xx, y0 = yy, x1 = xx2, y1 = yy, lwd = 10, col = cols)
text(rep(min(xx2), m), yy, paste(formatC(apply(DBR_area, 2 , max, na.rm = T),
                                         digits = 1, format = "f"), 
                 c("species", expression(km^2)), sep = " "), 
                 cex = 0.8, pos = 2)
dev.off()

########################

# analyses

pca <- prcomp(area_niche_DBR_fulldt_clean[, c("bio1", "bio2",
                                              "bio3", "bio4",
                                              "bio5", "bio6",
                                              "bio7", "bio8",
                                              "bio9", "bio10",
                                              "bio11", "bio12",
                                              "bio13", "bio14",
                                              "bio15", "bio16",
                                              "bio17", "bio18",
                                              "bio19", "aet2",
                                              "ai2", "alpha2",
                                              "etyr2", "npp08.18")], 
              center = TRUE, scale. = TRUE)

summary(pca, cutoff = 0.001, loadings = TRUE)

pcfactorloadings <- pca$rotation[, 1:3]

fulldt_analyses <- cbind(area_niche_DBR_fulldt_clean[, c("species", "area_sdm",
                                                         "speciesDBR", 
                                                         "genusDBR",
                                                         "familyDBR",
                                                         "orderDBR",
                                                         "classDBR",
                                                         "temp_niche_breadth",
                                                         "prec_niche_breadth")],
                        pca$x[, 1:3])

# verify if there is collinearity amog variables (OK)
result <- correlation(fulldt_analyses)
s <- summary(result)

pdf("figures/new/FigS1.pdf", width = 10)
plot(s)
dev.off()

# phylogenetic signal

fulldt_analyses2 <- fulldt_analyses %>% 
  dplyr::select(species, area_sdm, temp_niche_breadth, prec_niche_breadth,
                PC1, PC2, PC3) %>%
  remove_rownames() %>%
  column_to_rownames(var = "species")

p4d_tr <- phylo4d(tree, fulldt_analyses2)
phylosig_tr <- phyloSignal(p4d = p4d_tr, method = c("Lambda", "K"), 
                           reps = 1000)


# 1000 replicates
K <- K_p <- L <- L_p <- matrix(nrow = 6, ncol = 1000)
rownames(K) <- rownames(K_p) <- rownames(L) <- rownames(L_p) <- 
  colnames(fulldt_analyses2)
for (i in 1:1000) {
  p4d <- phylo4d(trees.pruned[[i]], fulldt_analyses2)
  phylosig <- phyloSignal(p4d = p4d, method = c("Lambda", "K"), reps = 1000)
  K[, i] <- phylosig$stat$K
  K_p[, i] <- phylosig$pvalue$K
  L[, i] <- phylosig$stat$Lambda
  L_p[, i] <- phylosig$pvalue$Lambda
}

res_phylosig <- data.frame(
    median_K = rowMeans(K),

    min_K = apply(K, 1, min),

    max_K = apply(K, 1, max),
 
    median_K_p = rowMeans(K_p),

    min_K_p = apply(K_p, 1, min),

    max_K_p = apply(K_p, 1, max),

    median_L = rowMeans(L),

    min_L = apply(L, 1, min),

    max_L = apply(L, 1, max),

    median_L_p = rowMeans(L_p),

    min_L_p = apply(L_p, 1, min),

    max_L_p = apply(L_p, 1, max)
)

res_phylosig2 <- res_phylosig %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  unite(K, c(min_K, max_K), sep = "-") %>%
  mutate(K = purrr::map(K, ~paste0("(", .x, ")"))) %>%
  unite(K_p, c(min_K_p, max_K_p), sep = "-") %>%
  mutate(K_p = purrr::map(K_p, ~paste0("(", .x, ")"))) %>%
  unite(L, c(min_L, max_L), sep = "-") %>%
  mutate(L = purrr::map(L, ~paste0("(", .x, ")"))) %>%
  unite(L_p, c(min_L_p, max_L_p), sep = "-") %>%
  mutate(L_p = purrr::map(L_p, ~paste0("(", .x, ")"))) %>%
  dplyr::select(K, K_p, L, L_p)

for (i in 1:nrow(res_phylosig2)) {
  res_phylosig2$K[i] <- paste0(round(phylosig_tr$stat$K[i], 3), 
                               " ", res_phylosig2$K[i])
  res_phylosig2$K_p[i] <- paste0(round(phylosig_tr$pvalue$K[i], 3), 
                               " ", res_phylosig2$K_p[i])
  res_phylosig2$L[i] <- paste0(round(phylosig_tr$stat$L[i], 3), 
                               " ", res_phylosig2$L[i])
  res_phylosig2$L_p[i] <- paste0(round(phylosig_tr$pvalue$L[i], 3), 
                               " ", res_phylosig2$L_p[i])
}

res_phylosig_final <- do.call(cbind, res_phylosig2)
rownames(res_phylosig_final) <- rownames(res_phylosig2)

write.csv(res_phylosig_final, "tables/new/Table1_unformatted.csv")

# pgls

layout(matrix(1:6, ncol = 3))
hist(fulldt_analyses$temp_niche_breadth)
hist(fulldt_analyses$prec_niche_breadth)
hist(fulldt_analyses$area)
hist(fulldt_analyses$PC1)
hist(fulldt_analyses$PC2)
hist(fulldt_analyses$PC3)

# scaled the predictors
fulldt_analyses$scale_nbtemp <- scale(fulldt_analyses$temp_niche_breadth, 
                                      center = TRUE)
fulldt_analyses$scale_nbprec <- scale(fulldt_analyses$prec_niche_breadth, 
                                      center = TRUE)
fulldt_analyses$scale_area <- scale(fulldt_analyses$area, center = TRUE)
fulldt_analyses$scale_pc1 <- scale(fulldt_analyses$PC1, center = TRUE)
fulldt_analyses$scale_pc2 <- scale(fulldt_analyses$PC2, center = TRUE)
fulldt_analyses$scale_pc3 <- scale(fulldt_analyses$PC3, center = TRUE)

layout(matrix(1:6, ncol = 3))
hist(fulldt_analyses$scale_nbtemp)
hist(fulldt_analyses$scale_nbprec)
hist(fulldt_analyses$scale_area)
hist(fulldt_analyses$scale_pc1)
hist(fulldt_analyses$scale_pc2)
hist(fulldt_analyses$scale_pc3)

comp_dat <- comparative.data(data = fulldt_analyses, phy = tree, 
                             names.col = "species", vcv.dim = 2, 
                             warn.dropped = TRUE)

comp_dat_trees <- list()
for (i in 1:length(trees)) {
  comp_dat_trees[[i]] <- comparative.data(data = fulldt_analyses, 
                                          phy = trees.pruned[[i]],
                                          names.col = "species", 
                                          vcv.dim = 2, 
                                          warn.dropped = TRUE)
}

# speciesDBR
res_pgls_spp <- pgls(speciesDBR ~ scale_area + scale_pc1 + scale_pc2 + 
                     scale_pc3 + scale_nbprec + scale_nbtemp, lambda = "ML",
                     data = comp_dat)

res_pgls_trees_spp <- list()
for (i in 1:length(comp_dat_trees)) {
  res_pgls_trees_spp[[i]] <- pgls(speciesDBR ~ scale_area + scale_pc1 +
                                  scale_pc2 + scale_pc3 + scale_nbprec + 
                                  scale_nbtemp, lambda = "ML",
                                  data = comp_dat_trees[[i]])
}

# genusDBR
res_pgls_gen <- pgls(genusDBR ~ scale_area + scale_pc1 + scale_pc2 + 
                     scale_pc3 + scale_nbprec + scale_nbtemp, lambda = "ML",
                     data = comp_dat)

res_pgls_trees_gen <- list()
for (i in 1:length(comp_dat_trees)) {
  res_pgls_trees_gen[[i]] <- pgls(genusDBR ~ scale_area + scale_pc1 +
                                  scale_pc2 + scale_pc3 + scale_nbprec + 
                                  scale_nbtemp, lambda = "ML",
                                  data = comp_dat_trees[[i]])
}

# familyDBR
res_pgls_fam <- pgls(familyDBR ~ scale_area + scale_pc1 + scale_pc2 + 
                     scale_pc3 + scale_nbprec + scale_nbtemp, lambda = "ML",
                     data = comp_dat)

res_pgls_trees_fam <- list()
for (i in 1:length(comp_dat_trees)) {
  res_pgls_trees_fam[[i]] <- pgls(familyDBR ~ scale_area + scale_pc1 +
                                  scale_pc2 + scale_pc3 + scale_nbprec + 
                                  scale_nbtemp, lambda = "ML",
                                  data = comp_dat_trees[[i]])
}

# orderDBR
res_pgls_ord <- pgls(orderDBR ~ scale_area + scale_pc1 + scale_pc2 + 
                     scale_pc3 + scale_nbprec + scale_nbtemp, lambda = "ML",
                     data = comp_dat)

res_pgls_trees_ord <- list()
for (i in 1:length(comp_dat_trees)) {
  res_pgls_trees_ord[[i]] <- pgls(orderDBR ~ scale_area + scale_pc1 +
                                  scale_pc2 + scale_pc3 + scale_nbprec + 
                                  scale_nbtemp, lambda = "ML",
                                  data = comp_dat_trees[[i]])
}

# classDBR
res_pgls_cla <- pgls(classDBR ~ scale_area + scale_pc1 + scale_pc2 + 
                     scale_pc3 + scale_nbprec + scale_nbtemp, lambda = "ML",
                     data = comp_dat)

res_pgls_trees_cla <- list()
for (i in 1:length(comp_dat_trees)) {
  res_pgls_trees_cla[[i]] <- pgls(classDBR ~ scale_area + scale_pc1 +
                                  scale_pc2 + scale_pc3 + scale_nbprec + 
                                  scale_nbtemp, lambda = "ML",
                                  data = comp_dat_trees[[i]])
}

stats_pgls <- function(data, data_trees) {
  
  final <- as.data.frame(matrix(ncol = 5, nrow = 6))
  for (i in 2:7) {
    
    for (j in 1:4) {
      
      stats <- numeric()
      for (k in 1:length(data_trees)) {
        stats[k] <- summary(data_trees[[k]])$coefficients[i, j]
      }

      range <- paste0("(", round(min(stats), 3), "-", round(max(stats), 3), ")")

      final[i-1, j] <- paste0(round(summary(data)$coefficients[i, j], 3), 
                              " ", range)

    }

    if (i == 7) {
      stats <- numeric()
      for (k in 1:length(data_trees)) {
        stats[k] <- summary(data_trees[[k]])$r.squared
      }
      
      range <- paste0("(", round(min(stats), 3), "-", round(max(stats), 3), ")")

      final[1, 5] <- paste0(round(summary(data)$r.squared, 3), 
                            " ", range)

    }

  }

  final <- cbind(rownames(summary(data)$coefficients)[2:7], final)
  colnames(final) <- c("Coefficient", colnames(summary(data)$coefficients),
                       "R_squared")

  return(final)
}

res_pgls <- data.frame(
  rbind(cbind(DBR = "speciesDBR", stats_pgls(res_pgls_spp, res_pgls_trees_spp)),
        cbind(DBR = "genusDBR", stats_pgls(res_pgls_gen, res_pgls_trees_gen)),
        cbind(DBR = "familyDBR", stats_pgls(res_pgls_fam, res_pgls_trees_fam)),
        cbind(DBR = "orderDBR", stats_pgls(res_pgls_ord, res_pgls_trees_ord)),
        cbind(DBR = "classDBR", stats_pgls(res_pgls_cla, res_pgls_trees_cla))))

write.csv(res_pgls, "tables/new/Table2_unformatted.csv")


