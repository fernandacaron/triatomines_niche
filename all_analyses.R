rm(list = ls())

setwd("~/Documents/lab/triatomines")

#library(rnaturalearth)
#library(rnaturalearthdata)
#library(ggspatial)
#library(sf)
#library(jsonlite)
#library(mapdata)
library(raster)
#library(mapview)
#library(ggplot2)
library(dplyr)
#library(maptools)
#library(rgdal)
#library(scrubr)
library(tidymodels)
library(tidyverse)
#library(purrr)
library(dismo)
#library(stringr)
library(rgbif)
library(tmap)
#library(hypervolume)
library(CoordinateCleaner)
library(stringi)

data <- read.csv("data/MasterTable11.csv", strip.white = TRUE, 
                 na.strings = c("", "NA"), stringsAsFactors = TRUE)

# species that cannot be used to calculate niche
subdat <- data.frame(KBGenus_Sp = c("Triatoma longipennis", 
                                    "Triatoma gusayana", 
                                    "Triatoma phyllosoma", 
                                    "Triatoma petrocchiae",
                                    "Psammolestes tertius",
                                    "Psammolestes coreodes", 
                                    "Rhodnius NA", 
                                    "Triatoma melanica"))

# dropping those factor levels
finaldt  <- data %>% anti_join(., subdat) %>% droplevels()

# function to fix names
fix_species_name <- function(column) {
  require(stringr)
  column_new <- stringr::str_replace(column, "(?s) .*", "sp")
  column_new
}

# host
finaldt$Host_Genus_sp <- str_replace(as.character(finaldt$Host_Genus_sp), "NA NA", 'NA')
finaldt$Host_Genus_sp <- str_replace(as.character(finaldt$Host_Genus_sp), "(?s) NA", " sp")
finaldt$Host_Genus_sp[which(finaldt$Host_Genus_sp == "NA rattus")] <- "Rattus rattus"
finaldt$Host_Genus_sp[which(finaldt$Host_Genus_sp == "Felis spp.")] <- "Felis sp"
finaldt$Host_Genus_sp <- as.factor(finaldt$Host_Genus_sp)
finaldt$Host_Genus_sp[which(finaldt$Host_Genus_sp == "NA")] <- NA
finaldt$Host_genus[which(is.na(finaldt$Host_genus) & !is.na(finaldt$Host_Genus_sp))] <- "Rattus"

# bugs
finaldt$KBGenus_Sp <- str_replace(as.character(finaldt$KBGenus_Sp), "NA NA", 'NA')
finaldt$KBGenus_Sp <- str_replace(as.character(finaldt$KBGenus_Sp), "(?s) NA", " sp")
finaldt$KBGenus_Sp <- as.factor(finaldt$KBGenus_Sp)

# country 
finaldt$Country <- str_replace(as.character(finaldt$Country), "Panamá", "Panama")
finaldt$Country <- str_replace(as.character(finaldt$Country), "Guyana Francesa", "French Guiana")
finaldt <- finaldt[!grepl("Trinidad and Tobago", finaldt$Country), ]

# estimating species-level and genus-level DBR for each species
# with NAs 
dietbreadthdt <- finaldt %>% 
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR = n_distinct(Family, na.rm = FALSE),
            orderDBR = n_distinct(Order, na.rm = FALSE),
            classDBR = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup()

# without NAs
dietbreadthdtnoNA <- finaldt %>% 
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
  na.omit() %>% filter(species != "Triatoma rubrofasciata")

dietbreadthdtnoNA_clean <- dietbreadthdtnoNA %>% ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR,
                genusDBR,
                familyDBR,
                orderDBR,
                classDBR) %>%
  na.omit() %>% filter(species != "Triatoma rubrofasciata")

genus_name_unique <- as.list(unique(sub('(^\\w+)\\s.+', '\\1', dietbreadthdt$KBGenus_Sp)))
species_name_unique <- as.list(unique(sub('(^\\w+)\\s(+)', '\\2', dietbreadthdt$KBGenus_Sp)))
species_name_unique <- unique(levels(dietbreadthdt$KBGenus_Sp))

# defininig the extent of latin america for cropping later
LatAmextent <- raster::extent(-125, -28, -70, 48) 
# same here, but with different name
LatAm.ext <- raster::extent(-125, -28, -70, 48)

gbif_taxon_keys <- as.character()
for (i in 1:length(species_name_unique)) {
  gbif_taxon_keys[i] <- name_backbone(species_name_unique[i])$usageKey
}

gbif_data <- occ_download(pred_in("taxonKey", gbif_taxon_keys),
                          pred("hasGeospatialIssue", FALSE),
                          pred("hasCoordinate", TRUE),
                          pred_not(pred_in("basisOfRecord", "FOSSIL_SPECIMEN")),
                          format = "SIMPLE_CSV")

gbif_download <- gbif_data %>%
  occ_download_get() %>%
  occ_download_import() %>%
  setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
  cc_sea() %>% # remove from ocean 
  cc_val() %>%
  cc_zero() %>%
  glimpse() # look at results of pipeline

# subsetting to only include species which we have the Diet Breadth
gbif_data_clean <- gbif_download[gbif_download$species %in% species_name_unique,]

# selecting the variables 
gbif_data_clean <- gbif_data_clean %>% dplyr::select(longitude = decimallongitude, latitude = decimallatitude, genus, species, country = countrycode)

# need to remove all NA's or function will not work
gbif_data_clean2 <- na.omit(gbif_data_clean) %>% droplevels(.)

write.csv(gbif_data_clean2, "data/Full_GBIF_data.csv")
#gbif_data_clean2 <- read.csv("data/Full_GBIF_data.csv")

# species into polygons
# now let's make our spatial data object just as we did in the last lesson
spatialdt <- SpatialPoints(coords = cbind(gbif_data_clean2$longitude, 
                                          gbif_data_clean2$latitude),
                           proj4string = CRS("+init=epsg:4326"))
spatialdt$species <- gbif_data_clean2$species
spatialdt$country <- gbif_data_clean2$country

# obtaining environmental data 
mlong <- mean(spatialdt@coords[,1])
mlat <- mean(spatialdt@coords[,2])
# query current environmental data for our location
clim.current <- getData("worldclim", var = "bio", res = 10, lat = mlat, 
                        lon = mlong)

# crop environmental variable extents
LatAm_clim <- crop(clim.current, LatAmextent)
LatAm_clim_df <- raster::extract(LatAm_clim, LatAmextent)

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
                     final_model = final_modsujarel,
                     model_performance = test_performance,
                     spatial_output_presence = sdm_spatial,
                     spatial_observed_data = spatial_presence)
  
  return(outputlist)

}

# 54.1124 mins
a <- Sys.time()
sdm_model_fulldt <- gbif_data_clean2 %>% 
  select(-country) %>% 
  nest(data = -species) %>%
  mutate(sdm = map(data, ~sdm_rf_manual(.x))) 
b <- Sys.time()
b-a
#There were issues with some computations   A: x1
#There were issues with some computations   A: x4
#There were issues with some computations   A: x2
#There were issues with some computations   A: x3
#There were issues with some computations   A: x2
#There were issues with some computations   A: x3
#There were issues with some computations   A: x1
#There were issues with some computations   A: x3
#There were issues with some computations   A: x1
#There were issues with some computations   A: x8

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
area_model_fulldt <- sdm_model_fulldt %>% select(-data) %>%  unnest(sdm) %>%
  group_by(species) %>% slice(4) %>%
  mutate(area_sdm = map(sdm, ~calculate_SDM_area(.x))) %>% select(-sdm) %>%
  unnest(area_sdm) %>% slice(3) %>% unnest(area_sdm)


## Getting the parameters from the model fit ##
param_random_forest_fulldt <- sdm_model_fulldt %>% select(-data) %>%  unnest(sdm) %>%
  group_by(species) %>% slice(3) %>% unnest() %>%
  select(-.estimator) %>%  tidyr::spread(., '.metric', '.estimate')

## Converting SpatialPointsDataFrame of the SDM output into a raster (for plotting) ##
raster_obj <- sdm_model_fulldt %>% select(-data) %>% unnest(sdm) %>% 
  group_by(species) %>% slice(4) %>% ungroup() %>% mutate(raster_obj = map(sdm, ~raster_obj(.x))) %>% 
  select(-sdm)


#There were issues with some computations   A: x1
#There were issues with some computations   A: x6
#There were issues with some computations   A: x1
#There were issues with some computations   A: x1
#There were issues with some computations   A: x3

calculate_SDM_area <- function(data, threshold = 0.5){
  unprojected_obj <- rasterFromXYZ(data)
  crs(unprojected_obj) <- projection(data)
  area_proj <- raster::area(unprojected_obj)
  avg_area_proj <- cellStats(unprojected_obj > threshold, sum) * area_proj
  output <- list(min_value = area_proj@data@min,
                 max_value = area_proj@data@max,
                 average_area_km2 = mean(getValues(avg_area_proj)))
  return(output)
}

## Estimating niche distribution based on SDM ##
area_model_fulldt <- sdm_model_fulldt %>% select(-data) %>%  unnest(sdm) %>%
  group_by(species) %>% slice(4) %>%
  mutate(area_sdm = map(sdm, ~calculate_SDM_area(.x))) %>% select(-sdm) %>%
  unnest(area_sdm) %>% slice(3) %>% unnest(area_sdm)

## Getting the parameters from the model fit ##
param_random_forest_fulldt <- sdm_model_fulldt %>% select(-data) %>%  unnest(sdm) %>%
  group_by(species) %>% slice(3) %>% unnest() %>%
  select(-.estimator) %>%  tidyr::spread(., '.metric', '.estimate')

### ~~~~~~~~ DIET BREADTH DATA ~~~~~~~~~~~ ####

## Calculating percentages of each host ##
dietbreadth_perc <- finaldt %>% dplyr::select(Study:Country, Lat:Host_genus, N_feeds) %>%
  group_by(species = KBGenus_Sp, host_family = Family, Method) %>% summarise(feeds = sum(N_feeds, na.rm = TRUE))

# making own palette
library("RColorBrewer")
my_cols <- c(brewer.pal(9, "Set1"), 
             brewer.pal(12, "Paired"), 
             brewer.pal(8, "Dark2"), 
             brewer.pal(8, "Accent"),
             brewer.pal(9, "Pastel1"),
             brewer.pal(8, "Pastel2"),
             brewer.pal(8, "Set2"),
             brewer.pal(12, "Set3"),
             brewer.pal(6, "RdGy"),
             brewer.pal(9, "BrBG"))


dietbreadth_perc %>% filter(!is.na(host_family)) %>%
  ggplot(aes(x = species, y = feeds, fill = host_family)) + ylab('Proportion of diet')+
  geom_bar(position = "fill", stat="identity") +coord_flip() + theme_bw() +
  guides(fill=FALSE) + #facet_grid(~Method) + 
  scale_fill_manual('Hosts', values = my_cols)# theme(axis.text.x = element_text(angle = 90)) 




##~~ join area with DBR ~~~
areaDBR_fulldt <- area_model_fulldt %>% inner_join(., dietbreadthdt_clean, by = 'species') #%>%
#separate(., species, into = c("Genus", "sp"), sep = " ")



### ~~~~~ NICHE BREADTH ~~~~~~~~ ####
## ~~ Predicted niche 
predicted_nichebreadth_fulldt <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>%  
  unnest(sdm) %>%
  group_by(species) %>% 
  slice(4) %>% 
  mutate(pred_nichebreadth = map(sdm, ~estimate_nichebreadth(.x))) %>%
  dplyr::select(-sdm) %>% unnest(pred_nichebreadth)


### ~~~>> putting all info into one dataset <<<~~~~ ####
area_niche_DBR_fulldt <- areaDBR_fulldt %>% 
  inner_join(., predicted_nichebreadth_fulldt, by = "species")












### ~~~~~~~~~ CALCULATING HYPERVOLUMES ~~~~~~~~~~ ####
# Estimating hypervolumes 
sdm_hypervol_estimates_fulldt <- sdm_model_fulldt %>% 
  dplyr::select(-data) %>%  
  unnest(sdm) %>%
  group_by(species) %>% 
  slice(4) %>%
  mutate(hypervol_calculation = map2(sdm, species, ~hypervolumefromSpatialPoints(.x, .y))) 

#Extracting information 
extrated_hyper_nested_fulldt <- sdm_hypervol_estimates_fulldt %>% 
  group_by(species) %>%
  mutate(hyper_var_imp = map(hypervol_calculation, 'Hypervolume_var_imp'),
         hyper_trimmed = map(hypervol_calculation, 'Hypervolume_trimmed')) %>%
  select(-hypervol_calculation) #%>% mutate(species1 = species) %>%
# separate(., species, into = c("Genus", "sp"), sep = " ")

## ~~~~~ Hypervolume volume ~~~~~~
extracted_hyper_volume <- extrated_hyper_nested_fulldt %>% 
  select(species, hyper_trimmed) %>%
  mutate(hyper_volume = map(hyper_trimmed, function(.x) {.x@Volume})) %>% 
  select(-hyper_trimmed) %>%
  unnest()



### ~~~~~>>>> All data combined <<<<<~~~~~~####
complete_DBRareanichehyper_dt <- inner_join(extracted_hyper_volume, area_niche_DBR_fulldt)
###~~~~~~~~ Writing up the info as new dataframes ~~~~~~~######
write.csv(complete_DBRareanichehyper_dt, "Full_Final_modelled_triatomidae_data_April2021.csv")
## ---------------- <<<>>>> --------------

## Getting the taxon Key ###
gbif_data_key <- data.frame()
for(i in 1:length(gbif_data)){
  if(i < 2){
    gbif_data_key <- data.frame(taxonKey = gbif_data[[i]]$taxonKey,
                                gbifID = gbif_data[[i]]$gbifID)
  } else {
    inter_key<- data.frame(taxonKey = gbif_data[[i]]$taxonKey,
                           gbifID = gbif_data[[i]]$gbifID)
    
    gbifdf <-bind_rows(gbif_data_key, inter_key)
  }
}

write.csv(gbif_data_key, "gbifKey_Triatomidae_data_April2021.csv")


































#########################################################################
######                                                              #####
###### ~~~~~~~~~~~~~~~~~~~ REGRESSION MODELS ~~~~~~~~~~~~~~~~~~~~~~ #####
######                                                              #####
#########################################################################
library(ape)
library(adephylo)
library(phylosignal)
library(ggrepel)
library(stringr)
library(nlme)

## Loading diet data (with NAs)
complete_DBRareanichehyper_dt <- read.csv("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Triatomidae-project/Full_Final_modelled_triatomidae_data_April2021.csv",
                                          header = TRUE)
# Removing the duplicated Rhodinius ecuadoriensis
complete_DBRareanichehyper_dt <- complete_DBRareanichehyper_dt[-18,]
## Counting missing values
nacount <- finaldt %>% 
  dplyr::select(KBGenus_Sp, Class:Host_genus) 
nacount_table <- data.frame(sum_NAs = apply(nacount, 2, function(x) sum(is.na(x))),
                            obs = nrow(finaldt))
nacount_table$prop_missing <- with(nacount_table, sum_NAs/obs * 100)
nacount_table







## Just assigning the previous full dataset to a new object for the sake of the script 
complete_DBRareanichehyper_dt_FULL <- complete_DBRareanichehyper_dt

### ~~~~~ Loading metadata with species name
metadata_tree <- read.csv("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Triatomidae-project/Completed_data/Meta_data_tree.csv",
                          strip.white = TRUE,
                          stringsAsFactors = FALSE)
## ~~~  Loading tree
tree <- ape::read.tree("/Users/s23jm9/Dropbox/University of Aberdeen/Research-Projects/1.Research Projects/Triatomidae-project/arbol_check.nex")

## Function to change labels based on metadata ###
change_tree_labels <- function(tree, metadata){
  new_tree <- c()
  tree_test <- tree
  metadata_tree <- metadata
  for(i in 1:length(tree_test$tip.label)){
    for(j in 1:nrow(metadata_tree)){
      tree_test$tip.label[i] <- ifelse(stringi::stri_compare(tree_test$tip.label[i], 
                                                             metadata_tree$Tip_label[j]) == 0,
                                       metadata_tree$names_full[j],
                                       next)
    }
  }
  new_tree <- tree_test
  return(new_tree)
}

## ~~~ Tree which is relabelled 
tree2 <- change_tree_labels(tree, metadata_tree)

## ~~~ Prunning phylo tree
labelstipstree <- data.frame(labels = tree2$tip.label)
labelstipsdf <- data.frame(labels = complete_DBRareanichehyper_dt_FULL$species)

testdata_tree <- inner_join(labelstipstree, labelstipsdf) %>% mutate(tip.label = labels) %>% filter(tip.label != "Trub") %>%
  #select(-labels) %>% 
  tibble::column_to_rownames(., var = "tip.label") 


nottree <- geiger::name.check(tree2, testdata_tree)$tree_not_data

### ~~~ >><<< Tree which is cleaned and relabelled ~~~ <<<>>> ###
tree_clean <-ape::drop.tip(tree2, nottree)
# ~~~~~~~~~<<<<<>>>>> ~~~~~~~~~ <<<>>>>> ~~~~~~~~~~~~ #

# Prunning the data to include only taxa in the tree
complete_DBRareanichehyper_dt_FULL_prunned <- complete_DBRareanichehyper_dt_FULL[complete_DBRareanichehyper_dt_FULL$species %in% tree_clean$tip.label,] %>%
  dplyr::mutate(species1 = species) %>%
  column_to_rownames("species1")
## ~~~~~~~~~~~~~~~~~~



## Correlation plot DBR metrics ###
corrplot::corrplot(as.matrix(cor(complete_DBRareanichehyper_dt_FULL[6:9])), 
                   method= c("number"),
                   type = "lower",
                   tl.col="black",
                   tl.srt=60)

### Correlation plot clim vars ###
corrplot::corrplot(as.matrix(cor(na.omit(LatAm_clim_df[,c("bio1", "bio4", "bio7", "bio12", "bio15")]))), 
                   method= c("number"),
                   type = "lower",
                   tl.col="black",
                   tl.srt=60)




##### ~~~ Regression Models ~~~ #####
### ~~~~~~~~~ Full Model  ~~~~~~~~~~~ ######
## Correlations
with(complete_DBRareanichehyper_dt_FULL_prunned,
     cor.test(hyper_volume, pred_nichebreadth))
with(complete_DBRareanichehyper_dt_FULL_prunned,
     cor.test(area_sdm, pred_nichebreadth))
with(complete_DBRareanichehyper_dt_FULL_prunned,
     cor.test(area_sdm, hyper_volume))

with(complete_DBRareanichehyper_dt_FULL,
     cor.test(hyper_volume, pred_nichebreadth))
# ~~ accounting for collinearity between pred_niche and hypervol
complete_DBRareanichehyper_dt_FULL_prunned$res_hyper <- residuals(lm(complete_DBRareanichehyper_dt_FULL_prunned$hyper_volume ~ complete_DBRareanichehyper_dt_FULL_prunned$pred_nichebreadth))
complete_DBRareanichehyper_dt_FULL$res_hyper <- residuals(lm(complete_DBRareanichehyper_dt_FULL$hyper_volume ~ complete_DBRareanichehyper_dt_FULL$pred_nichebreadth))


#~ Phylogenetically controlled
fullmodel_area0 <-nlme::gls(log(familyDBR) ~ log(area_sdm) +
                              pred_nichebreadth + res_hyper, 
                             data=complete_DBRareanichehyper_dt_FULL_prunned, 
                             correlation = ape::corPagel(1, tree_clean))

summary(fullmodel_area0)
anova(fullmodel_area0)



#~ LM
fullmodel_LM <- lm(log(familyDBR) ~ log(area_sdm) +
                              pred_nichebreadth + res_hyper, 
                            data=complete_DBRareanichehyper_dt_FULL)
summary(fullmodel_LM)
anova(fullmodel_LM)







#### ~~~~~ PLOTTING  ~~~~~~ ###########
## ~~ (Phylogenetically controlled) ~~~~
pred_fullmodel <- data.frame(complete_DBRareanichehyper_dt_FULL_prunned, 
                        phylo_cont = exp(predict(fullmodel_area0, complete_DBRareanichehyper_dt_FULL_prunned)))
predframe_fullmodel  <- with(complete_DBRareanichehyper_dt_FULL_prunned,data.frame(familyDBR,
                                                                              wow=pred_fullmodel$phylo_cont,
                                                                              lwr=exp(log(pred_fullmodel$phylo_cont)-1.96*0.08),
                                                                              upr=exp(log(pred_fullmodel$phylo_cont)+1.96*0.08)))


fullplotdt_fullmodel<- data.frame(area_sdm = complete_DBRareanichehyper_dt_FULL_prunned$area_sdm,
                                  niche_fam = complete_DBRareanichehyper_dt_FULL_prunned$pred_nichebreadth,
                                  res_hyper_fam = complete_DBRareanichehyper_dt_FULL_prunned$res_hyper,
                              DBR_fam = complete_DBRareanichehyper_dt_FULL_prunned$familyDBR,
                              fitline = predframe_fullmodel$wow,
                              lwr = predframe_fullmodel$lwr,
                              upr = predframe_fullmodel$upr,
                              species = complete_DBRareanichehyper_dt_FULL_prunned$species) %>% 
  mutate(species1 = species) %>%
  separate(., species1, into = c("Genus", "sp"), sep = " ")


# plot 
phylotree_plot_DBR_geo0 <- ggplot(fullplotdt_fullmodel) + 
  geom_point(aes(x = log(area_sdm), 
                 y =  log(DBR_fam),
                 fill = Genus),
             pch = 21,
             size = 3, alpha = 0.5) +
  geom_smooth(aes(y = log(fitline), x = log(area_sdm)), size = 1, method = "lm", se = FALSE, col = "black") +
  theme_bw() + ggtitle('Phylogenetically controlled') +
  theme(plot.title = element_text(face = "bold.italic", size = 9),
        panel.grid = element_blank(),
        legend.position = "none") +
  xlab(expression('log(Predicted area) in '~km^2)) +
  ylab('log(Family-level DBR)') +
  scale_fill_manual(breaks = c("Mepraia", "Panstrongylus", "Rhodnius", "Triatoma", "Paratriatoma"),
                    values = c("orangered3", "royalblue4", "green3", "steelblue2", "brown2"))

phylotree_plot_DBR_geo0






phylotree_plot_DBR_niche <- ggplot(fullplotdt_fullmodel ) + 
  geom_point(aes(x = niche_fam, 
                 y =  log(DBR_fam),
                 fill = Genus),
             pch = 21,
             size = 3, alpha = 0.5) +
  geom_smooth(aes(y = log(fitline), x = niche_fam), size = 1, method = "lm", se = FALSE, col = "black") +
  theme_bw() + ggtitle('Phylogenetically controlled') +
  theme(plot.title = element_text(face = "bold.italic", size = 9),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab(expression('log(Family-level DBR')) +
  xlab('Niche breadth') +
  scale_fill_manual(breaks = c("Mepraia", "Panstrongylus", "Rhodnius", "Triatoma", "Paratriatoma"),
                    values = c("orangered3", "royalblue4", "green3", "steelblue2", "brown2"))

phylotree_plot_DBR_niche





## ~~ LM
## ~~ Plotting (LM) ~~~~
pred_fullmodelLM <- data.frame(complete_DBRareanichehyper_dt_FULL, 
                             LM_cont = exp(predict(fullmodel_LM, complete_DBRareanichehyper_dt_FULL)))
predframe_fullmodel_LM  <- with(complete_DBRareanichehyper_dt_FULL,data.frame(familyDBR,
                                                                                   wow=pred_fullmodelLM$LM_cont))


fullplotdt_fullmodel_LM <- data.frame(area_sdm = complete_DBRareanichehyper_dt_FULL$area_sdm,
                                  niche_fam = complete_DBRareanichehyper_dt_FULL$pred_nichebreadth,
                                  res_hyper_fam = complete_DBRareanichehyper_dt_FULL$res_hyper,
                                  DBR_fam = complete_DBRareanichehyper_dt_FULL$familyDBR,
                                  fitline = predframe_fullmodel_LM$wow,
                                  species = complete_DBRareanichehyper_dt_FULL$species) %>% 
  mutate(species1 = species) %>%
  separate(., species1, into = c("Genus", "sp"), sep = " ")



# plot
notree_plot_DBR_geo2 <- ggplot(fullplotdt_fullmodel_LM ) + 
  geom_jitter(aes(x = log(DBR_fam), 
                  y = log(area_sdm),
                  fill = Genus), width = 0.2,
              pch = 21,
              size = 3, alpha = 0.5) +
  geom_smooth(aes(x = log(DBR_fam), 
                  y = log(area_sdm)), 
              method = "lm", se = FALSE, col = "black") +
  theme_bw() + ggtitle('Linear regression') +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold.italic", size = 9),
        panel.grid = element_blank()) +
  xlab(expression('log(Predicted area) in '~km^2)) +
  ylab('log(Family-level DBR)') + 
  scale_fill_manual(breaks = c("Mepraia", "Panstrongylus", "Rhodnius", "Triatoma",
                               "Belminus", "Cavernicola", "Meccus", "Paratriatoma",
                               "Psammolestes", "Eratyrus"),
                    values = c("orangered3", "royalblue4", "green3", "steelblue2",
                               "darkorange2", "mediumorchid2", "goldenrod1", "brown2", 
                               "darkgreen", "red2"))
notree_plot_DBR_geo2






notree_plot_DBR_niche <- ggplot(fullplotdt_fullmodel_LM ) + 
  geom_jitter(aes(y = log(DBR_fam), 
                  x = niche_fam,
                  fill = Genus), width = 0.2,
              pch = 21,
              size = 3, alpha = 0.5) +
  geom_smooth(aes(y = log(DBR_fam), 
                  x = niche_fam), 
              method = "lm", se = FALSE, col = "black") +
  theme_bw() + ggtitle('Linear regression') +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold.italic", size = 9),
        panel.grid = element_blank()) +
  xlab(expression('Niche breadth')) +
  ylab('log(Family-level DBR)') + 
  scale_fill_manual(breaks = c("Mepraia", "Panstrongylus", "Rhodnius", "Triatoma",
                               "Belminus", "Cavernicola", "Meccus", "Paratriatoma",
                               "Psammolestes", "Eratyrus"),
                    values = c("orangered3", "royalblue4", "green3", "steelblue2",
                               "darkorange2", "mediumorchid2", "goldenrod1", "brown2", 
                               "darkgreen", "red2"))
notree_plot_DBR_niche


## Multipanel
library(patchwork)

DBRplot <- (phylotree_plot_DBR_geo0 | notree_plot_DBR_geo2) / (phylotree_plot_DBR_niche | notree_plot_DBR_niche )
DBRplot <- DBRplot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14))

DBRplot # Saved as APRIL2021_FINAL_AreavsDBR.pdf














##### ~~~ Feeding trajectories  ~~~~~~~
# data frame
Diet_trajectoriesdt <- complete_DBRareanichehyper_dt_FULL %>%
  dplyr::select(species, speciesDBR:classDBR) %>%
  pivot_longer(cols = speciesDBR:classDBR,
               names_to = "Level",
               values_to = "DBR") %>%
  arrange(desc(Level, DBR)) %>%
  tidyr::separate(col = "species", into = c("Genus", "Sp"), sep = " ", remove = FALSE)


# Converting first letter to capital and removing the DBR suffix
Diet_trajectoriesdt$Level <- as.factor(stringr::str_to_title(stringr::str_replace_all(Diet_trajectoriesdt$Level,
                                                                                      "DBR", "")))

## Reordering for plotting 
Diet_trajectoriesdt$Level <- relevel(Diet_trajectoriesdt$Level, ref = "Class")
Diet_trajectoriesdt$Level <- relevel(Diet_trajectoriesdt$Level, ref = "Order")
Diet_trajectoriesdt$Level <- relevel(Diet_trajectoriesdt$Level, ref = "Family")
Diet_trajectoriesdt$Level <- relevel(Diet_trajectoriesdt$Level, ref = "Genus")
Diet_trajectoriesdt$Level <- relevel(Diet_trajectoriesdt$Level, ref = "Species")
 

# Calculating mean trajectories per Genus
Mean_diettrajectory_genus <- Diet_trajectoriesdt %>%
  group_by(Genus, Level) %>% 
  summarise(DBR_genus = mean(DBR))

# Calculating overall mean trajectory
Mean_diettrajectory_all <- Diet_trajectoriesdt %>%
  group_by(Level) %>% 
  summarise(DBR_all = mean(DBR))

ggplot(Diet_trajectoriesdt) + 
  geom_line(aes(x = Level, y = DBR, group = species, colour = Genus), alpha = 0.3, size = 0.4, col = "grey60") +
  geom_line(data = Mean_diettrajectory_genus, aes(x = Level, y = DBR_genus, group = Genus, col= Genus), size = 0.8, alpha = 0.6) +
  #geom_line(data = Mean_diettrajectory_all, aes(x = Level, y = DBR_all, group = 1), size = 1.2, col = "black", alpha = 0.6, lty = "dashed") +
  scale_x_discrete(expand = c(0.01,0.03)) +
  theme_linedraw() + 
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.justification = "left",
        legend.title = element_text(family = "Arial", size = 12, face = "bold")) +
  xlab(expression('Host taxonomic level')) +
  ylab('Number of hosts') + 
  scale_colour_manual(breaks = c("Mepraia", "Panstrongylus", "Rhodnius", "Triatoma",
                               "Belminus", "Cavernicola", "Meccus", "Paratriatoma",
                               "Psammolestes", "Eratyrus"),
                    values = c("orangered3", "royalblue4", "green3", "steelblue2",
                               "darkorange2", "mediumorchid2", "goldenrod1", "brown2", 
                               "darkgreen", "red2"))
  














### ~~~~ Diet analyses ~~~~~~ ######
library("RColorBrewer")
my_cols <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(3, "Accent"))

DBR_full <-  finaldt %>% 
  dplyr::select(Study:Country, Lat:Host_genus, N_feeds) %>%
  dplyr::filter(N_feeds > 0) %>%
  dplyr::arrange(desc(N_feeds)) 

DBR_full$Class <- as.factor(DBR_full$Class)
DBR_full$Class <- ifelse(DBR_full$Class == "Sauropsida (Reptilia)", "Reptilia",  as.character(DBR_full$Class))
table(DBR_full$Class)




### Pareto function
pareto.MLE <- function(data, i) { 
  require(dplyr)
  X <- data[i]
  n <- length(X)
  m <- min(X)
  a <- n/sum(log(X)-log(m))
  return( c(m,a) ) 
}

library(boot)
pareto.MLE(complete_DBRareanichehyper_dt_FULL$familyDBR)
pareto_boot_sp <- boot(complete_DBRareanichehyper_dt_FULL$familyDBR, pareto.MLE, R = 1000)
boot.ci(pareto_boot_sp, index = 2)






library("RColorBrewer")
my_cols <- c(brewer.pal(8, "Set2"), 
             brewer.pal(12, "Paired"), 
             brewer.pal(8, "Dark2"), 
             brewer.pal(3, "Accent"),
             brewer.pal(8, "Set1"))

DBR_full  %>% 
  ggplot(aes(x = KBGenus_Sp, y = N_feeds, fill = Class)) + ylab('Proportion of feeds on host')+ xlab('Species') +
  geom_bar(position = "fill", stat="identity") +coord_flip() + theme_bw() + #+guides(fill=FALSE)
  scale_fill_manual('Hosts', values = my_cols) 
#+ coord_polar("y", start=0) # theme(axis.text.x = element_text(angle = 90)) 


### Table with % of feed of each class ####
DBR_full$Class <- as.factor(DBR_full$Class)
percentage_classDBR <- DBR_full %>% dplyr::select(KBGenus_Sp, Class, N_feeds) %>% 
  group_by(KBGenus_Sp, Class) %>% summarise(N_feeds2 = sum(N_feeds))  %>%
  group_by(KBGenus_Sp) %>% mutate(totalfeeds = round(N_feeds2/sum(N_feeds2)*100, 2)) %>% 
  dplyr::select(-N_feeds2) %>%
  pivot_wider(KBGenus_Sp, names_from = Class, values_from = totalfeeds, 
              values_fill = 0)





DBR_full$Family <- as.factor(DBR_full$Family)
DBR_full$Family <- forcats::fct_explicit_na(DBR_full$Family)

percentage_famDBR <- DBR_full %>% dplyr::select(KBGenus_Sp, Family, N_feeds) %>% 
  group_by(KBGenus_Sp, Family) %>% summarise(N_feeds2 = sum(N_feeds))  %>%
  group_by(KBGenus_Sp) %>% mutate(totalfeeds = round(N_feeds2/sum(N_feeds2)*100, 2)) %>% 
  dplyr::select(-N_feeds2) %>%
  pivot_wider(KBGenus_Sp, names_from = Family, values_from = totalfeeds, 
              values_fill = 0)

percentage_famDBR

write.csv(percentage_classDBR, "Table1a_MainText.csv") ##Writing the table##
write.csv(percentage_famDBR, "Table1b_MainText.csv")






### ~~~ Analysing human blood feeding vs niche ~~~ ####
percHuman <- data.frame(species = percentage_famDBR$KBGenus_Sp,
                        Human_perc = percentage_famDBR$`Human Homininae`)

humandt <- inner_join(complete_DBRareanichehyper_dt_FULL, percHuman)

ggplot(humandt, aes(x = genusDBR, y = Human_perc)) + 
  geom_point() +
  geom_smooth(se = FALSE, method = "gam")
