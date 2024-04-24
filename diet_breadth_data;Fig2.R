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

# estimating species-level and genus-level DBR for each species
dietbreadthdt_clean <- data %>% 
  dplyr::select(KBTribe:Country, Lat:Host_genus, N_feeds) %>%
  group_by(KBTribe, KBGenus_Sp, KBGenus) %>% 
  summarise(speciesDBR = n_distinct(Host_Genus_sp, na.rm = FALSE),
            genusDBR = n_distinct(Host_genus, na.rm = FALSE),
            familyDBR = n_distinct(Family, na.rm = FALSE),
            orderDBR = n_distinct(Order, na.rm = FALSE),
            classDBR = n_distinct(Class, na.rm = FALSE)) %>%
  ungroup() %>% 
  dplyr::select(species = KBGenus_Sp,
                speciesDBR,
                genusDBR,
                familyDBR,
                orderDBR,
                classDBR) %>%
  na.omit()

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

write.csv(dietbreadthdt_clean, "data/dietbreadthdt_clean.csv")
