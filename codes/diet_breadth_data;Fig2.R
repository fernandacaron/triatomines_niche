rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(stringr)
library(dplyr)

data <- read.csv("data/MasterTable11.csv", strip.white = TRUE, 
                 na.strings = c("", "NA"), stringsAsFactors = TRUE)

data$KBTribe[data$KBGenus_Sp == "Rhodnius ecuadoriensis"] <- "Rhodniini"

data$Host_Genus_sp <- str_replace(as.character(data$Host_Genus_sp), "NA NA", 'NA')
data$Host_Genus_sp <- str_replace(as.character(data$Host_Genus_sp), "(?s) NA", " sp")
data$Host_Genus_sp[which(data$Host_Genus_sp == "NA rattus")] <- "Rattus rattus"
data$Host_Genus_sp[which(data$Host_Genus_sp == "Felis spp.")] <- "Felis sp"
data$Host_Genus_sp <- as.factor(data$Host_Genus_sp)
data$Host_Genus_sp[which(data$Host_Genus_sp == "NA")] <- NA
data$Host_genus[which(is.na(data$Host_genus) & !is.na(data$Host_Genus_sp))] <- "Rattus"

data$KBGenus_Sp <- str_replace(as.character(data$KBGenus_Sp), "NA NA", 'NA')
data$KBGenus_Sp <- str_replace(as.character(data$KBGenus_Sp), "(?s) NA", " sp")
data$KBGenus_Sp <- as.factor(data$KBGenus_Sp)

data$Country <- str_replace(as.character(data$Country), "Panamá", "Panama")
data$Country <- str_replace(as.character(data$Country), "Guyana Francesa", "French Guiana")

data <- cbind(data, NA, NA, NA)
colnames(data)[29:31] <- c("Domiciliary", "Peridomiciliary", "Sylvatic")
for (i in 1:nrow(data)) {
  spp_i <- data$KBGenus_Sp[i]
  
  habitat_i <- as.character(data$Habitat[i])
  
  list_habitat <- unique(do.call(c, lapply(habitat_i, 
                                           function(x) strsplit(x, "_")[[1]])))
  
  if ("Domiciliary" %in% list_habitat) {
    data$Domiciliary[i] <- 1
  } else {
    data$Domiciliary[i] <- 0
  }
  
  if ("Peridomiciliary" %in% list_habitat) {
    data$Peridomiciliary[i] <- 1
  } else {
    data$Peridomiciliary[i] <- 0
  }
  
  if ("Sylvatic" %in% list_habitat) {
    data$Sylvatic[i] <- 1
  } else {
    data$Sylvatic[i] <- 0
  }
  
}

# estimating species-level and genus-level DBR for each species
dietbreadthdt_clean <- data %>% 
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR_all = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR_all = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR_all = n_distinct(Family, na.rm = FALSE),
            orderDBR_all = n_distinct(Order, na.rm = FALSE),
            classDBR_all = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR_all,
                genusDBR_all,
                familyDBR_all,
                orderDBR_all,
                classDBR_all) %>%
  na.omit()

dietbreadthdt_dom_clean <- data %>%
  filter(Domiciliary == 1) %>%
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR_dom = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR_dom = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR_dom = n_distinct(Family, na.rm = FALSE),
            orderDBR_dom = n_distinct(Order, na.rm = FALSE),
            classDBR_dom = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR_dom,
                genusDBR_dom,
                familyDBR_dom,
                orderDBR_dom,
                classDBR_dom) %>%
  na.omit()

dietbreadthdt_per_clean <- data %>%
  filter(Peridomiciliary == 1) %>%
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR_per = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR_per = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR_per = n_distinct(Family, na.rm = FALSE),
            orderDBR_per = n_distinct(Order, na.rm = FALSE),
            classDBR_per = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR_per,
                genusDBR_per,
                familyDBR_per,
                orderDBR_per,
                classDBR_per) %>%
  na.omit()

dietbreadthdt_syl_clean <- data %>%
  filter(Sylvatic == 1) %>%
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR_syl = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR_syl = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR_syl = n_distinct(Family, na.rm = FALSE),
            orderDBR_syl = n_distinct(Order, na.rm = FALSE),
            classDBR_syl = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR_syl,
                genusDBR_syl,
                familyDBR_syl,
                orderDBR_syl,
                classDBR_syl) %>%
  na.omit()

habitat <- data$Habitat
names(habitat) <- data$KBGenus_Sp
dietbreadthdt_all_clean <- dietbreadthdt_clean
dietbreadthdt_all_clean <- cbind(dietbreadthdt_all_clean, 
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
colnames(dietbreadthdt_all_clean)[7:24] <- c(colnames(dietbreadthdt_dom_clean)[2:6],
                                             colnames(dietbreadthdt_per_clean)[2:6],
                                             colnames(dietbreadthdt_syl_clean)[2:6],
                                             "Domiciliary",
                                             "Peridomiciliary",
                                             "Sylvatic")

for (i in 1:nrow(dietbreadthdt_all_clean)) {
  spp_i <- dietbreadthdt_all_clean$species[i]
  
  if (spp_i %in% dietbreadthdt_dom_clean$species) {
    dietbreadthdt_all_clean$speciesDBR_dom[i] <- 
      dietbreadthdt_dom_clean$speciesDBR_dom[dietbreadthdt_dom_clean$species == spp_i]
    
    dietbreadthdt_all_clean$genusDBR_dom[i] <- 
      dietbreadthdt_dom_clean$genusDBR_dom[dietbreadthdt_dom_clean$species == spp_i]
    
    dietbreadthdt_all_clean$familyDBR_dom[i] <- 
      dietbreadthdt_dom_clean$familyDBR_dom[dietbreadthdt_dom_clean$species == spp_i]
    
    dietbreadthdt_all_clean$orderDBR_dom[i] <- 
      dietbreadthdt_dom_clean$orderDBR_dom[dietbreadthdt_dom_clean$species == spp_i]
    
    dietbreadthdt_all_clean$classDBR_dom[i] <- 
      dietbreadthdt_dom_clean$classDBR_dom[dietbreadthdt_dom_clean$species == spp_i]
  } else {
    dietbreadthdt_all_clean$speciesDBR_dom[i] <- NA
    dietbreadthdt_all_clean$genusDBR_dom[i] <- NA
    dietbreadthdt_all_clean$familyDBR_dom[i] <- NA
    dietbreadthdt_all_clean$orderDBR_dom[i] <- NA
    dietbreadthdt_all_clean$classDBR_dom[i] <- NA
  }
  
  if (spp_i %in% dietbreadthdt_per_clean$species) {
    dietbreadthdt_all_clean$speciesDBR_per[i] <- 
      dietbreadthdt_per_clean$speciesDBR_per[dietbreadthdt_per_clean$species == spp_i]
    
    dietbreadthdt_all_clean$genusDBR_per[i] <- 
      dietbreadthdt_per_clean$genusDBR_per[dietbreadthdt_per_clean$species == spp_i]
    
    dietbreadthdt_all_clean$familyDBR_per[i] <- 
      dietbreadthdt_per_clean$familyDBR_per[dietbreadthdt_per_clean$species == spp_i]
    
    dietbreadthdt_all_clean$orderDBR_per[i] <- 
      dietbreadthdt_per_clean$orderDBR_per[dietbreadthdt_per_clean$species == spp_i]
    
    dietbreadthdt_all_clean$classDBR_per[i] <- 
      dietbreadthdt_per_clean$classDBR_per[dietbreadthdt_per_clean$species == spp_i]
  } else {
    dietbreadthdt_all_clean$speciesDBR_per[i] <- NA
    dietbreadthdt_all_clean$genusDBR_per[i] <- NA
    dietbreadthdt_all_clean$familyDBR_per[i] <- NA
    dietbreadthdt_all_clean$orderDBR_per[i] <- NA
    dietbreadthdt_all_clean$classDBR_per[i] <- NA
  }
  
  if (spp_i %in% dietbreadthdt_syl_clean$species) {
    dietbreadthdt_all_clean$speciesDBR_syl[i] <- 
      dietbreadthdt_syl_clean$speciesDBR_syl[dietbreadthdt_syl_clean$species == spp_i]
    
    dietbreadthdt_all_clean$genusDBR_syl[i] <- 
      dietbreadthdt_syl_clean$genusDBR_syl[dietbreadthdt_syl_clean$species == spp_i]
    
    dietbreadthdt_all_clean$familyDBR_syl[i] <- 
      dietbreadthdt_syl_clean$familyDBR_syl[dietbreadthdt_syl_clean$species == spp_i]
    
    dietbreadthdt_all_clean$orderDBR_syl[i] <- 
      dietbreadthdt_syl_clean$orderDBR_syl[dietbreadthdt_syl_clean$species == spp_i]
    
    dietbreadthdt_all_clean$classDBR_syl[i] <- 
      dietbreadthdt_syl_clean$classDBR_syl[dietbreadthdt_syl_clean$species == spp_i]
  } else {
    dietbreadthdt_all_clean$speciesDBR_syl[i] <- NA
    dietbreadthdt_all_clean$genusDBR_syl[i] <- NA
    dietbreadthdt_all_clean$familyDBR_syl[i] <- NA
    dietbreadthdt_all_clean$orderDBR_syl[i] <- NA
    dietbreadthdt_all_clean$classDBR_syl[i] <- NA
  }
 
  habitat_i <- as.character(unique(habitat[names(habitat) == spp_i]))
  
  list_habitat <- unique(do.call(c, lapply(habitat_i, 
                                           function(x) strsplit(x, "_")[[1]])))
  
  if ("Domiciliary" %in% list_habitat) {
    dietbreadthdt_all_clean$Domiciliary[dietbreadthdt_all_clean$species == spp_i] <- 1
  } else {
    dietbreadthdt_all_clean$Domiciliary[dietbreadthdt_all_clean$species == spp_i] <- 0
  }
  
  if ("Peridomiciliary" %in% list_habitat) {
    dietbreadthdt_all_clean$Peridomiciliary[dietbreadthdt_all_clean$species == spp_i] <- 1
  } else {
    dietbreadthdt_all_clean$Peridomiciliary[dietbreadthdt_all_clean$species == spp_i] <- 0
  }
  
  if ("Sylvatic" %in% list_habitat) {
    dietbreadthdt_all_clean$Sylvatic[dietbreadthdt_all_clean$species == spp_i] <- 1
  } else {
    dietbreadthdt_all_clean$Sylvatic[dietbreadthdt_all_clean$species == spp_i] <- 0
  }
  
}

# calculating percentages of each host
dietbreadth_perc <- data %>%
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(species = KBGenus_Sp, host_family = Family) %>% 
  summarise(feeds = sum(N_feeds, na.rm = TRUE)) %>% 
  droplevels()

# Figure 2
pdf("figures/Figure2_unformatted.pdf", height = 10)
dietbreadth_perc %>% 
  filter(!is.na(host_family)) %>%
  ggplot(aes(x = species, y = feeds, fill = host_family)) + 
  ylab("Proportion of diet") +
  geom_bar(position = "fill", stat = "identity") + 
  coord_flip() + 
  theme_bw(base_line_size = 0, base_rect_size = 0) +
  guides(fill = FALSE) + 
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "none") + 
  geom_text(aes(y = feeds, label = toupper(as.character(host_family))), 
            size = 2, position = position_fill(vjust = 0.5), color = "white")
dev.off()

write.csv(dietbreadthdt_all_clean, "data/dietbreadthdt_all_clean.csv")
