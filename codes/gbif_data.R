rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(dplyr)
library(rgbif)
library(CoordinateCleaner)

dietbreadthdt_clean <- read.csv("data/dietbreadthdt_clean.csv", row.names = 1)

genus_name_unique <- as.list(unique(sub('(^\\w+)\\s.+', '\\1', dietbreadthdt_clean$species)))
species_name_unique <- as.list(unique(sub('(^\\w+)\\s(+)', '\\2', dietbreadthdt_clean$species)))
species_name_unique <- unique(dietbreadthdt_clean$species)

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

# Table S1
write.csv(gbif_data_clean2, "data/Full_GBIF_data.csv", row.names = F)
