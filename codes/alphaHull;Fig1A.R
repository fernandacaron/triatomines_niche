rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(geodata)
library(raster)
library(dplyr)
library(tidyr)
library(rangeBuilder)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(phytools)
library(stringi)

gbif_data <- read.csv("data/Full_GBIF_data.csv")
gbif_data$genus <- stri_split_fixed(gbif_data$species, " ", simplify = T)[, 1]

gbif_data_clean <- gbif_data

gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Eratyrus cuspidatus" &
                                     gbif_data_clean$country == "BR"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Mepraia spinolai" &
                                     gbif_data_clean$country == "MX"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Panstrongylus chinai" &
                                     gbif_data_clean$country %in% c("BR", "VE")), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Panstrongylus geniculatus" &
                                       gbif_data_clean$country == "US"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Panstrongylus megistus" &
                                       gbif_data_clean$country == "MX"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Rhodnius neglectus" &
                                       gbif_data_clean$country == "MX"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Rhodnius prolixus" &
                                       gbif_data_clean$country %in% c("US", "BR", "GF")), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma dimidiata" &
                                       gbif_data_clean$country == "BR"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma infestans" &
                                       gbif_data_clean$country %in% c("US", "VE", "MX")), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma rubida" &
                                       gbif_data_clean$longitude == -81.37944), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma rubrofasciata" &
                                       gbif_data_clean$country == "US"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma rubrovaria" &
                                       gbif_data_clean$country == "VE"), ]
gbif_data_clean <- gbif_data_clean[!(gbif_data_clean$species == "Triatoma vitticeps" &
                                       gbif_data_clean$country %in% c("VE", "MX")), ]

#write.csv(gbif_data_clean, "tables/occ.csv", row.names = F)

# Table S2
doi <- read.csv("tables/doi_gbif.csv")$doi

tab_s1 <- gbif_data_clean
spp <- unique(gbif_data_clean$species)
for (i in 1:length(spp)) {
  tab_s1$doi[tab_s1$species == spp[i]] <- doi[i]
}

write.csv(tab_s2, "tables/table_s2.csv", row.names = F)

tr <- read.nexus("data/tr_cal_f.tre")
gbif_data_clean2 <- gbif_data_clean[gbif_data_clean$species %in% stri_replace_all_fixed(tr$tip.label, "_", " "), ]

gbif_data_clean2$genus[gbif_data_clean2$species == "Meccus pallidipennis"] <- "Triatoma"
gbif_data_clean2$species[gbif_data_clean2$species == "Meccus pallidipennis"] <- "Triatoma pallidipennis"

gbif_data_clean2$genus[gbif_data_clean2$species == "Triatoma lecticularia"] <- "Paratriatoma"
gbif_data_clean2$species[gbif_data_clean2$species == "Triatoma lecticularia"] <- "Paratriatoma lecticularia"

gbif_data_clean2$genus[gbif_data_clean2$species == "Triatoma tibiamaculata"] <- "Panstrongylus"
gbif_data_clean2$species[gbif_data_clean2$species == "Triatoma tibiamaculata"] <- "Panstrongylus tibiamaculatus"

gbif_sf <- sf::st_as_sf(gbif_data_clean2, coords = c("longitude", "latitude"),
                        crs = "+proj=longlat +datum=WGS84")

map <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf", 
                    continent = c("south america", "north america", "central america"))

colors <- magma(12)[c(1, 5, 7, 8, 9, 10, 11, 12)]
names(colors) <- c("Triatoma", "Panstrongylus", "Rhodnius", 
                   "Psammolestes", "Mepraia", "Paratriatoma", "Cavernicola")

pdf("figures/spatial_data.pdf", width = 10, height = 7)
ggplot() +
  geom_sf(data = map, color = "gray", fill = "lightgray") +
  xlim(-180, 0) + 
  geom_sf(data = gbif_sf %>% filter(genus == "Triatoma"), 
          aes(color = "Triatoma")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Panstrongylus"), 
          aes(color = "Panstrongylus")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Rhodnius"), 
          aes(color = "Rhodnius")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Psammolestes"), 
          aes(color = "Psammolestes")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Mepraia"), 
          aes(color = "Mepraia")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Paratriatoma"), 
          aes(color = "Paratriatoma")) +
  geom_sf(data = gbif_sf %>% filter(genus == "Cavernicola"), 
          aes(color = "Cavernicola")) +
  scale_color_manual(values = colors) + 
  theme_void()  
dev.off()

trgetAlphaHull <- function(data) {
  dups <- duplicated(round(data[c("longitude", "latitude")], 2))
  data <- data[!dups, ]
  ahull <- NULL
  try({ahull <- getDynamicAlphaHull(data,
                               fraction = 1, partCount = 1, buff = 10000,
                               initialAlpha = 0.1, 
                               coordHeaders = c("longitude", "latitude"),
                               clipToCoast = "terrestrial",
                               alphaIncrement = 1, verbose = F)
  })
  return(ahull)
}

alp_hull <- gbif_data_clean %>% 
  dplyr::select(-country) %>% 
  nest(data = -species) %>%
  mutate(coord = purrr::map(data, ~dplyr::select(.x, -genus))) %>%
  mutate(ahull = list(NULL))

# manually removing problematic points
spp_i <- alp_hull$coord[[59]]
ahull_i <- getAlphaHull(spp_i)
#spp_j <- spp_i[c(1:767, 769:840, 842:843, 845:955, 957:1309, 1311:4585),]
#ahull_i <- getAlphaHull(spp_j)
#save(alp_hull, file = "data/alp_hull.RData")

#load(file = "data/alp_hull.RData")
alp_hull <- alp_hull %>% 
  mutate(area = purrr::map(ahull, function(.x) (sum(st_area(.x[[1]]))) / 1000000))

# query current environmental data for our location
clim_current <- worldclim_global("bio", res = 10, path = tempdir())
names(clim_current) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                         "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
                         "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
clim_current <- stack(clim_current)

# raster with the same dimensions as above
r <- raster(ncols = 2160, nrows = 1080, ymn = -90, vals = 1:(2160*1080))

extract_env <- function(data, clim = clim_current) {
  
  values_env <- raster::extract(clim_current, data)
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
  
  temp_niche_breadth <- (max(values_env$bio5, na.rm = T)) - 
    (min(values_env$bio6, na.rm = T))
  names(temp_niche_breadth) <- "temp_niche_breadth"
  prec_niche_breadth <- max(values_env$bio16, na.rm = T) - 
    min(values_env$bio17, na.rm = T)
  names(prec_niche_breadth) <- "prec_niche_breadth"
  
  return(c(env_avg, temp_niche_breadth, prec_niche_breadth))
}

clim_data <- alp_hull %>% 
  mutate(env_vals = purrr::map(coord, ~extract_env(.x))) 

dietbreadth <- read.csv("data/dietbreadthdt_all_clean.csv", row.names = 1)

full_data <- clim_data %>% 
  dplyr::select(c(species, area, env_vals)) %>%
  unnest(area) %>%
  unnest_wider(env_vals) %>%
  inner_join(., dietbreadth, by = "species") 

write.csv(full_data, "data/dietbreadth_niche_habitat_11apr24.csv", 
          row.names = F)
