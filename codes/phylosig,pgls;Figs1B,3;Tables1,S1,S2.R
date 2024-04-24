rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(stringi)
library(phytools)
library(geiger)
library(patchwork)
library(phylosignal)
library(phylobase)
library(caper)
library(ggtree)
library(deeptime)
library(tibble)
library(ggplot2)

dat <- read.csv("data/dietbreadth_niche_habitat_11apr24.csv")

# loading phylogenies (MCC and 1000 random trees from posterior distribution)
tr <- read.nexus("data/tr_cal_f.tre")
tr_post <- read.nexus("data/1k_random_tria.trees")

# eliminate the samples not in the phylogeny and vice-versa
dat$species <- stri_replace_all_fixed(dat$species, " ", "_")
dat <- dat[dat$species %in% tr$tip.label, ]

tr <- drop.tip(tr, setdiff(tr$tip.label, dat$species))
trees <- list()
for (i in 1:1000) {
  trees[[i]] <- drop.tip(tr_post[[i]], setdiff(tr_post[[i]]$tip.label, dat$species))
}

pdf("figures/Figure1_scale.pdf")
revts(ggtree(tr_ladd) + geom_tiplab(size = 3)) +
  coord_geo(xlim = c(-60, 30), ylim = c(-2, 55), neg = TRUE, abbrv = F, dat = "epochs") +
  scale_x_continuous(breaks = seq(-80, 0, 20), labels = abs(seq(-80, 0, 20))) +
  theme_tree2()
dev.off()

dbr <- dat$speciesDBR
names(dbr) <- dat$species
dbr <- dbr[tr$tip.label]
dbr <- dbr/max(dbr)

area <- dat$area
names(area) <- dat$species
area <- area[tr$tip.label]
area <- area/max(area)

temp <- dat$temp_niche_breadth
names(temp) <- dat$species
temp <- temp[tr$tip.label]
temp <- temp/max(temp)

prec <- dat$prec_niche_breadth
names(prec) <- dat$species
prec <- prec[tr$tip.label]
prec <- prec/max(prec)

tr_ladd <- ladderize(tr, right = FALSE)
tr_ladd <- ips::fixNodes(tr_ladd)

cols <- c("#36877a", "#3a7c89", "#235d72","#123f5a")

pdf("figures/Figure1_tree.pdf", width = 10)
layout(matrix(1:5, ncol = 5), widths = c(5/10, 2/10, 2/10, 2/10, 2/10))
plotTree.barplot(tr_ladd, area, add = TRUE, width = 2, 
                 args.barplot = list(xlab = "Range area", col = cols[1], 
                                     border = cols[1], mar = c(5.1, 0, 2.1, 2.1),
                                     xlim = c(0, 1.4)))
plotTree.barplot(tr_ladd, dbr, add = TRUE, width = 2, 
                 args.barplot = list(xlab = "Diet niche breadth", col = cols[2], 
                                     border = cols[2], mar = c(5.1, 0 ,2.1, 2.1), 
                                     xlim = c(0, 1.4)), 
                 args.plotTree = list(plot = FALSE))
plotTree.barplot(tr_ladd, temp, add = TRUE,
                 args.barplot = list(xlab = "Temperature niche breadth", col = cols[3], 
                                     border = cols[3], mar = c(5.1, 0, 2.1, 2.1),
                                     xlim = c(0, 1.4)), 
                 args.plotTree = list(plot = FALSE))
plotTree.barplot(tr_ladd, prec, add = TRUE,
                 args.barplot = list(xlab = "Precipitation niche breadth", col = cols[4], 
                                     border = cols[4], mar = c(5.1, 0, 2.1, 2.1),
                                     xlim = c(0, 1.4)), 
                 args.plotTree = list(plot = FALSE))

dev.off()

# phylogenetic signal

dat_physig <- dat %>% 
  dplyr::select(species, area, temp_niche_breadth, prec_niche_breadth,
                speciesDBR, genusDBR, familyDBR, orderDBR, classDBR) %>%
  remove_rownames() %>%
  column_to_rownames(var = "species")

K <- K_p <- L <- L_p <- matrix(nrow = 1000, ncol = 8)
colnames(K) <- colnames(K_p) <- colnames(L) <- colnames(L_p) <- 
  colnames(dat_physig)
for (i in 1:1000) {
  p4d <- phylo4d(trees[[i]], dat_physig)
  phylosig <- phyloSignal(p4d = p4d, method = c("Lambda", "K"), reps = 1000)
  K[i, ] <- phylosig$stat$K
  K_p[i, ] <- phylosig$pvalue$K
  L[i, ] <- phylosig$stat$Lambda
  L_p[i, ] <- phylosig$pvalue$Lambda
}

res_physig <- cbind(
  paste0(apply(K, 2, function(x) round(median(x), 3)), " (", 
         apply(K, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(K, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(K_p, 2, function(x) round(median(x), 3)), " (", 
         apply(K_p, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(K_p, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(L, 2, function(x) round(median(x), 3)), " (", 
         apply(L, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(L, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(L_p, 2, function(x) round(median(x), 3)), " (", 
         apply(L_p, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(L_p, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")")
)

res_physig <- as.data.frame(cbind(colnames(K), res_physig))
colnames(res_physig) <- c("trait", "K", "K_p", "L", "L_p")

write.csv(res_physig, "tables/Table_physig.csv")

# simple regreesion
reg_spp_all <- lm(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat)
reg_gen_all <- lm(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat)
reg_fam_all <- lm(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat)
reg_ord_all <- lm(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat)
reg_cla_all <- lm(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat)

reg_spp_dom <- lm(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Domiciliary == 1, ])
reg_gen_dom <- lm(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Domiciliary == 1, ])
reg_fam_dom <- lm(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Domiciliary == 1, ])
reg_ord_dom <- lm(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Domiciliary == 1, ])
reg_cla_dom <- lm(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Domiciliary == 1, ])

reg_spp_per <- lm(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Peridomiciliary == 1, ])
reg_gen_per <- lm(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Peridomiciliary == 1, ])
reg_fam_per <- lm(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Peridomiciliary == 1, ])
reg_ord_per <- lm(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Peridomiciliary == 1, ])
reg_cla_per <- lm(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Peridomiciliary == 1, ])

reg_spp_syl <- lm(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Sylvatic == 1, ])
reg_gen_syl <- lm(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Sylvatic == 1, ])
reg_fam_syl <- lm(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Sylvatic == 1, ])
reg_ord_syl <- lm(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Sylvatic == 1, ])
reg_cla_syl <- lm(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) + 
                    sqrt(prec_niche_breadth), data = dat[dat$Sylvatic == 1, ])

res_reg <- data.frame(
  habitat = rep(c("All", "Domiciliary", "Peridomiciliary", "Sylvatic"), each = 5),
  
  DBR = rep(c("species", "genus", "family", "order", "class"), 4),
  
  Intercept = c(round(summary(reg_spp_all)$coefficients[1, 1], 3),
                round(summary(reg_gen_all)$coefficients[1, 1], 3),
                round(summary(reg_fam_all)$coefficients[1, 1], 3),
                round(summary(reg_ord_all)$coefficients[1, 1], 3),
                round(summary(reg_cla_all)$coefficients[1, 1], 3),
                round(summary(reg_spp_dom)$coefficients[1, 1], 3),
                round(summary(reg_gen_dom)$coefficients[1, 1], 3),
                round(summary(reg_fam_dom)$coefficients[1, 1], 3),
                round(summary(reg_ord_dom)$coefficients[1, 1], 3),
                round(summary(reg_cla_dom)$coefficients[1, 1], 3),
                round(summary(reg_spp_per)$coefficients[1, 1], 3),
                round(summary(reg_gen_per)$coefficients[1, 1], 3),
                round(summary(reg_fam_per)$coefficients[1, 1], 3),
                round(summary(reg_ord_per)$coefficients[1, 1], 3),
                round(summary(reg_cla_per)$coefficients[1, 1], 3),
                round(summary(reg_spp_syl)$coefficients[1, 1], 3),
                round(summary(reg_gen_syl)$coefficients[1, 1], 3),
                round(summary(reg_fam_syl)$coefficients[1, 1], 3),
                round(summary(reg_ord_syl)$coefficients[1, 1], 3),
                round(summary(reg_cla_syl)$coefficients[1, 1], 3)),
  
  DBR_estimate = c(round(summary(reg_spp_all)$coefficients[2, 1], 3),
                   round(summary(reg_gen_all)$coefficients[2, 1], 3),
                   round(summary(reg_fam_all)$coefficients[2, 1], 3),
                   round(summary(reg_ord_all)$coefficients[2, 1], 3),
                   round(summary(reg_cla_all)$coefficients[2, 1], 3),
                   round(summary(reg_spp_dom)$coefficients[2, 1], 3),
                   round(summary(reg_gen_dom)$coefficients[2, 1], 3),
                   round(summary(reg_fam_dom)$coefficients[2, 1], 3),
                   round(summary(reg_ord_dom)$coefficients[2, 1], 3),
                   round(summary(reg_cla_dom)$coefficients[2, 1], 3),
                   round(summary(reg_spp_per)$coefficients[2, 1], 3),
                   round(summary(reg_gen_per)$coefficients[2, 1], 3),
                   round(summary(reg_fam_per)$coefficients[2, 1], 3),
                   round(summary(reg_ord_per)$coefficients[2, 1], 3),
                   round(summary(reg_cla_per)$coefficients[2, 1], 3),
                   round(summary(reg_spp_syl)$coefficients[2, 1], 3),
                   round(summary(reg_gen_syl)$coefficients[2, 1], 3),
                   round(summary(reg_fam_syl)$coefficients[2, 1], 3),
                   round(summary(reg_ord_syl)$coefficients[2, 1], 3),
                   round(summary(reg_cla_syl)$coefficients[2, 1], 3)),
  
  DBR_Std_Error = c(round(summary(reg_spp_all)$coefficients[2, 2], 3),
                    round(summary(reg_gen_all)$coefficients[2, 2], 3),
                    round(summary(reg_fam_all)$coefficients[2, 2], 3),
                    round(summary(reg_ord_all)$coefficients[2, 2], 3),
                    round(summary(reg_cla_all)$coefficients[2, 2], 3),
                    round(summary(reg_spp_dom)$coefficients[2, 2], 3),
                    round(summary(reg_gen_dom)$coefficients[2, 2], 3),
                    round(summary(reg_fam_dom)$coefficients[2, 2], 3),
                    round(summary(reg_ord_dom)$coefficients[2, 2], 3),
                    round(summary(reg_cla_dom)$coefficients[2, 2], 3),
                    round(summary(reg_spp_per)$coefficients[2, 2], 3),
                    round(summary(reg_gen_per)$coefficients[2, 2], 3),
                    round(summary(reg_fam_per)$coefficients[2, 2], 3),
                    round(summary(reg_ord_per)$coefficients[2, 2], 3),
                    round(summary(reg_cla_per)$coefficients[2, 2], 3),
                    round(summary(reg_spp_syl)$coefficients[2, 2], 3),
                    round(summary(reg_gen_syl)$coefficients[2, 2], 3),
                    round(summary(reg_fam_syl)$coefficients[2, 2], 3),
                    round(summary(reg_ord_syl)$coefficients[2, 2], 3),
                    round(summary(reg_cla_syl)$coefficients[2, 2], 3)),
  
  DBR_t = c(round(summary(reg_spp_all)$coefficients[2, 3], 3),
            round(summary(reg_gen_all)$coefficients[2, 3], 3),
            round(summary(reg_fam_all)$coefficients[2, 3], 3),
            round(summary(reg_ord_all)$coefficients[2, 3], 3),
            round(summary(reg_cla_all)$coefficients[2, 3], 3),
            round(summary(reg_spp_dom)$coefficients[2, 3], 3),
            round(summary(reg_gen_dom)$coefficients[2, 3], 3),
            round(summary(reg_fam_dom)$coefficients[2, 3], 3),
            round(summary(reg_ord_dom)$coefficients[2, 3], 3),
            round(summary(reg_cla_dom)$coefficients[2, 3], 3),
            round(summary(reg_spp_per)$coefficients[2, 3], 3),
            round(summary(reg_gen_per)$coefficients[2, 3], 3),
            round(summary(reg_fam_per)$coefficients[2, 3], 3),
            round(summary(reg_ord_per)$coefficients[2, 3], 3),
            round(summary(reg_cla_per)$coefficients[2, 3], 3),
            round(summary(reg_spp_syl)$coefficients[2, 3], 3),
            round(summary(reg_gen_syl)$coefficients[2, 3], 3),
            round(summary(reg_fam_syl)$coefficients[2, 3], 3),
            round(summary(reg_ord_syl)$coefficients[2, 3], 3),
            round(summary(reg_cla_syl)$coefficients[2, 3], 3)),
  
  DBR_p = c(round(summary(reg_spp_all)$coefficients[2, 4], 3),
            round(summary(reg_gen_all)$coefficients[2, 4], 3),
            round(summary(reg_fam_all)$coefficients[2, 4], 3),
            round(summary(reg_ord_all)$coefficients[2, 4], 3),
            round(summary(reg_cla_all)$coefficients[2, 4], 3),
            round(summary(reg_spp_dom)$coefficients[2, 4], 3),
            round(summary(reg_gen_dom)$coefficients[2, 4], 3),
            round(summary(reg_fam_dom)$coefficients[2, 4], 3),
            round(summary(reg_ord_dom)$coefficients[2, 4], 3),
            round(summary(reg_cla_dom)$coefficients[2, 4], 3),
            round(summary(reg_spp_per)$coefficients[2, 4], 3),
            round(summary(reg_gen_per)$coefficients[2, 4], 3),
            round(summary(reg_fam_per)$coefficients[2, 4], 3),
            round(summary(reg_ord_per)$coefficients[2, 4], 3),
            round(summary(reg_cla_per)$coefficients[2, 4], 3),
            round(summary(reg_spp_syl)$coefficients[2, 4], 3),
            round(summary(reg_gen_syl)$coefficients[2, 4], 3),
            round(summary(reg_fam_syl)$coefficients[2, 4], 3),
            round(summary(reg_ord_syl)$coefficients[2, 4], 3),
            round(summary(reg_cla_syl)$coefficients[2, 4], 3)),
  
  Temp_estimate = c(round(summary(reg_spp_all)$coefficients[3, 1], 3),
                    round(summary(reg_gen_all)$coefficients[3, 1], 3),
                    round(summary(reg_fam_all)$coefficients[3, 1], 3),
                    round(summary(reg_ord_all)$coefficients[3, 1], 3),
                    round(summary(reg_cla_all)$coefficients[3, 1], 3),
                    round(summary(reg_spp_dom)$coefficients[3, 1], 3),
                    round(summary(reg_gen_dom)$coefficients[3, 1], 3),
                    round(summary(reg_fam_dom)$coefficients[3, 1], 3),
                    round(summary(reg_ord_dom)$coefficients[3, 1], 3),
                    round(summary(reg_cla_dom)$coefficients[3, 1], 3),
                    round(summary(reg_spp_per)$coefficients[3, 1], 3),
                    round(summary(reg_gen_per)$coefficients[3, 1], 3),
                    round(summary(reg_fam_per)$coefficients[3, 1], 3),
                    round(summary(reg_ord_per)$coefficients[3, 1], 3),
                    round(summary(reg_cla_per)$coefficients[3, 1], 3),
                    round(summary(reg_spp_syl)$coefficients[3, 1], 3),
                    round(summary(reg_gen_syl)$coefficients[3, 1], 3),
                    round(summary(reg_fam_syl)$coefficients[3, 1], 3),
                    round(summary(reg_ord_syl)$coefficients[3, 1], 3),
                    round(summary(reg_cla_syl)$coefficients[3, 1], 3)),
  
  Temp_Std_Error = c(round(summary(reg_spp_all)$coefficients[3, 2], 3),
                     round(summary(reg_gen_all)$coefficients[3, 2], 3),
                     round(summary(reg_fam_all)$coefficients[3, 2], 3),
                     round(summary(reg_ord_all)$coefficients[3, 2], 3),
                     round(summary(reg_cla_all)$coefficients[3, 2], 3),
                     round(summary(reg_spp_dom)$coefficients[3, 2], 3),
                     round(summary(reg_gen_dom)$coefficients[3, 2], 3),
                     round(summary(reg_fam_dom)$coefficients[3, 2], 3),
                     round(summary(reg_ord_dom)$coefficients[3, 2], 3),
                     round(summary(reg_cla_dom)$coefficients[3, 2], 3),
                     round(summary(reg_spp_per)$coefficients[3, 2], 3),
                     round(summary(reg_gen_per)$coefficients[3, 2], 3),
                     round(summary(reg_fam_per)$coefficients[3, 2], 3),
                     round(summary(reg_ord_per)$coefficients[3, 2], 3),
                     round(summary(reg_cla_per)$coefficients[3, 2], 3),
                     round(summary(reg_spp_syl)$coefficients[3, 2], 3),
                     round(summary(reg_gen_syl)$coefficients[3, 2], 3),
                     round(summary(reg_fam_syl)$coefficients[3, 2], 3),
                     round(summary(reg_ord_syl)$coefficients[3, 2], 3),
                     round(summary(reg_cla_syl)$coefficients[3, 2], 3)),
  
  Temp_t = c(round(summary(reg_spp_all)$coefficients[3, 3], 3),
             round(summary(reg_gen_all)$coefficients[3, 3], 3),
             round(summary(reg_fam_all)$coefficients[3, 3], 3),
             round(summary(reg_ord_all)$coefficients[3, 3], 3),
             round(summary(reg_cla_all)$coefficients[3, 3], 3),
             round(summary(reg_spp_dom)$coefficients[3, 3], 3),
             round(summary(reg_gen_dom)$coefficients[3, 3], 3),
             round(summary(reg_fam_dom)$coefficients[3, 3], 3),
             round(summary(reg_ord_dom)$coefficients[3, 3], 3),
             round(summary(reg_cla_dom)$coefficients[3, 3], 3),
             round(summary(reg_spp_per)$coefficients[3, 3], 3),
             round(summary(reg_gen_per)$coefficients[3, 3], 3),
             round(summary(reg_fam_per)$coefficients[3, 3], 3),
             round(summary(reg_ord_per)$coefficients[3, 3], 3),
             round(summary(reg_cla_per)$coefficients[3, 3], 3),
             round(summary(reg_spp_syl)$coefficients[3, 3], 3),
             round(summary(reg_gen_syl)$coefficients[3, 3], 3),
             round(summary(reg_fam_syl)$coefficients[3, 3], 3),
             round(summary(reg_ord_syl)$coefficients[3, 3], 3),
             round(summary(reg_cla_syl)$coefficients[3, 3], 3)),
  
  Temp_p = c(round(summary(reg_spp_all)$coefficients[3, 4], 3),
             round(summary(reg_gen_all)$coefficients[3, 4], 3),
             round(summary(reg_fam_all)$coefficients[3, 4], 3),
             round(summary(reg_ord_all)$coefficients[3, 4], 3),
             round(summary(reg_cla_all)$coefficients[3, 4], 3),
             round(summary(reg_spp_dom)$coefficients[3, 4], 3),
             round(summary(reg_gen_dom)$coefficients[3, 4], 3),
             round(summary(reg_fam_dom)$coefficients[3, 4], 3),
             round(summary(reg_ord_dom)$coefficients[3, 4], 3),
             round(summary(reg_cla_dom)$coefficients[3, 4], 3),
             round(summary(reg_spp_per)$coefficients[3, 4], 3),
             round(summary(reg_gen_per)$coefficients[3, 4], 3),
             round(summary(reg_fam_per)$coefficients[3, 4], 3),
             round(summary(reg_ord_per)$coefficients[3, 4], 3),
             round(summary(reg_cla_per)$coefficients[3, 4], 3),
             round(summary(reg_spp_syl)$coefficients[3, 4], 3),
             round(summary(reg_gen_syl)$coefficients[3, 4], 3),
             round(summary(reg_fam_syl)$coefficients[3, 4], 3),
             round(summary(reg_ord_syl)$coefficients[3, 4], 3),
             round(summary(reg_cla_syl)$coefficients[3, 4], 3)),
  
  Prec_estimate = c(round(summary(reg_spp_all)$coefficients[4, 1], 3),
                    round(summary(reg_gen_all)$coefficients[4, 1], 3),
                    round(summary(reg_fam_all)$coefficients[4, 1], 3),
                    round(summary(reg_ord_all)$coefficients[4, 1], 3),
                    round(summary(reg_cla_all)$coefficients[4, 1], 3),
                    round(summary(reg_spp_dom)$coefficients[4, 1], 3),
                    round(summary(reg_gen_dom)$coefficients[4, 1], 3),
                    round(summary(reg_fam_dom)$coefficients[4, 1], 3),
                    round(summary(reg_ord_dom)$coefficients[4, 1], 3),
                    round(summary(reg_cla_dom)$coefficients[4, 1], 3),
                    round(summary(reg_spp_per)$coefficients[4, 1], 3),
                    round(summary(reg_gen_per)$coefficients[4, 1], 3),
                    round(summary(reg_fam_per)$coefficients[4, 1], 3),
                    round(summary(reg_ord_per)$coefficients[4, 1], 3),
                    round(summary(reg_cla_per)$coefficients[4, 1], 3),
                    round(summary(reg_spp_syl)$coefficients[4, 1], 3),
                    round(summary(reg_gen_syl)$coefficients[4, 1], 3),
                    round(summary(reg_fam_syl)$coefficients[4, 1], 3),
                    round(summary(reg_ord_syl)$coefficients[4, 1], 3),
                    round(summary(reg_cla_syl)$coefficients[4, 1], 3)),
  
  Prec_Std_Error = c(round(summary(reg_spp_all)$coefficients[4, 2], 3),
                     round(summary(reg_gen_all)$coefficients[4, 2], 3),
                     round(summary(reg_fam_all)$coefficients[4, 2], 3),
                     round(summary(reg_ord_all)$coefficients[4, 2], 3),
                     round(summary(reg_cla_all)$coefficients[4, 2], 3),
                     round(summary(reg_spp_dom)$coefficients[4, 2], 3),
                     round(summary(reg_gen_dom)$coefficients[4, 2], 3),
                     round(summary(reg_fam_dom)$coefficients[4, 2], 3),
                     round(summary(reg_ord_dom)$coefficients[4, 2], 3),
                     round(summary(reg_cla_dom)$coefficients[4, 2], 3),
                     round(summary(reg_spp_per)$coefficients[4, 2], 3),
                     round(summary(reg_gen_per)$coefficients[4, 2], 3),
                     round(summary(reg_fam_per)$coefficients[4, 2], 3),
                     round(summary(reg_ord_per)$coefficients[4, 2], 3),
                     round(summary(reg_cla_per)$coefficients[4, 2], 3),
                     round(summary(reg_spp_syl)$coefficients[4, 2], 3),
                     round(summary(reg_gen_syl)$coefficients[4, 2], 3),
                     round(summary(reg_fam_syl)$coefficients[4, 2], 3),
                     round(summary(reg_ord_syl)$coefficients[4, 2], 3),
                     round(summary(reg_cla_syl)$coefficients[4, 2], 3)),
  
  Prec_t = c(round(summary(reg_spp_all)$coefficients[4, 3], 3),
             round(summary(reg_gen_all)$coefficients[4, 3], 3),
             round(summary(reg_fam_all)$coefficients[4, 3], 3),
             round(summary(reg_ord_all)$coefficients[4, 3], 3),
             round(summary(reg_cla_all)$coefficients[4, 3], 3),
             round(summary(reg_spp_dom)$coefficients[4, 3], 3),
             round(summary(reg_gen_dom)$coefficients[4, 3], 3),
             round(summary(reg_fam_dom)$coefficients[4, 3], 3),
             round(summary(reg_ord_dom)$coefficients[4, 3], 3),
             round(summary(reg_cla_dom)$coefficients[4, 3], 3),
             round(summary(reg_spp_per)$coefficients[4, 3], 3),
             round(summary(reg_gen_per)$coefficients[4, 3], 3),
             round(summary(reg_fam_per)$coefficients[4, 3], 3),
             round(summary(reg_ord_per)$coefficients[4, 3], 3),
             round(summary(reg_cla_per)$coefficients[4, 3], 3),
             round(summary(reg_spp_syl)$coefficients[4, 3], 3),
             round(summary(reg_gen_syl)$coefficients[4, 3], 3),
             round(summary(reg_fam_syl)$coefficients[4, 3], 3),
             round(summary(reg_ord_syl)$coefficients[4, 3], 3),
             round(summary(reg_cla_syl)$coefficients[4, 3], 3)),
  
  Prec_p = c(round(summary(reg_spp_all)$coefficients[4, 4], 3),
             round(summary(reg_gen_all)$coefficients[4, 4], 3),
             round(summary(reg_fam_all)$coefficients[4, 4], 3),
             round(summary(reg_ord_all)$coefficients[4, 4], 3),
             round(summary(reg_cla_all)$coefficients[4, 4], 3),
             round(summary(reg_spp_dom)$coefficients[4, 4], 3),
             round(summary(reg_gen_dom)$coefficients[4, 4], 3),
             round(summary(reg_fam_dom)$coefficients[4, 4], 3),
             round(summary(reg_ord_dom)$coefficients[4, 4], 3),
             round(summary(reg_cla_dom)$coefficients[4, 4], 3),
             round(summary(reg_spp_per)$coefficients[4, 4], 3),
             round(summary(reg_gen_per)$coefficients[4, 4], 3),
             round(summary(reg_fam_per)$coefficients[4, 4], 3),
             round(summary(reg_ord_per)$coefficients[4, 4], 3),
             round(summary(reg_cla_per)$coefficients[4, 4], 3),
             round(summary(reg_spp_syl)$coefficients[4, 4], 3),
             round(summary(reg_gen_syl)$coefficients[4, 4], 3),
             round(summary(reg_fam_syl)$coefficients[4, 4], 3),
             round(summary(reg_ord_syl)$coefficients[4, 4], 3),
             round(summary(reg_cla_syl)$coefficients[4, 4], 3)),
  
  R_squared = c(round(summary(reg_spp_all)$r.squared, 3),
                round(summary(reg_gen_all)$r.squared, 3),
                round(summary(reg_fam_all)$r.squared, 3),
                round(summary(reg_ord_all)$r.squared, 3),
                round(summary(reg_cla_all)$r.squared, 3),
                round(summary(reg_spp_dom)$r.squared, 3),
                round(summary(reg_gen_dom)$r.squared, 3),
                round(summary(reg_fam_dom)$r.squared, 3),
                round(summary(reg_ord_dom)$r.squared, 3),
                round(summary(reg_cla_dom)$r.squared, 3),
                round(summary(reg_spp_per)$r.squared, 3),
                round(summary(reg_gen_per)$r.squared, 3),
                round(summary(reg_fam_per)$r.squared, 3),
                round(summary(reg_ord_per)$r.squared, 3),
                round(summary(reg_cla_per)$r.squared, 3),
                round(summary(reg_spp_syl)$r.squared, 3),
                round(summary(reg_gen_syl)$r.squared, 3),
                round(summary(reg_fam_syl)$r.squared, 3),
                round(summary(reg_ord_syl)$r.squared, 3),
                round(summary(reg_cla_syl)$r.squared, 3))
)

write.csv(res_reg, "tables/Table_reg.csv", row.names = F)

p1 <- ggplot() +
  geom_point(data = dat, aes(x = log(speciesDBR), y = log(area)), size = 2) + 
  geom_smooth(data = dat, aes(x = log(speciesDBR), y = log(area)), method = "lm", 
              se = FALSE, linewidth = 1, linetype = 2, colour = "black") +
  geom_smooth(data = dat[dat$Domiciliary == 1, ], aes(x = log(speciesDBR), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 2, colour = "steelblue3") +
  geom_smooth(data = dat[dat$Peridomiciliary == 1, ], aes(x = log(speciesDBR), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 2, colour = "orange2") +
  geom_smooth(data = dat[dat$Sylvatic == 1, ], aes(x = log(speciesDBR), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 2, colour = "chartreuse3") +
  labs(x = "Diet niche breadth", y = expression(paste("log Range size (km"^"2", ")"))) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  guides(colour = FALSE) 

p2 <- ggplot() +
  geom_point(data = dat, aes(x = sqrt(temp_niche_breadth), y = log(area)), size = 2) + 
  geom_smooth(data = dat, aes(x = sqrt(temp_niche_breadth), y = log(area)), method = "lm", 
              se = FALSE, linewidth = 1, linetype = 1, colour = "black") +
  geom_smooth(data = dat[dat$Domiciliary == 1, ], aes(x = sqrt(temp_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 1, colour = "steelblue3") +
  geom_smooth(data = dat[dat$Peridomiciliary == 1, ], aes(x = sqrt(temp_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 1, colour = "orange2") +
  geom_smooth(data = dat[dat$Sylvatic == 1, ], aes(x = sqrt(temp_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 1, colour = "chartreuse3") +
  labs(x = "Temperature niche breadth", y = expression(paste("log Range size (km"^"2", ")"))) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  guides(colour = FALSE)

p3 <- ggplot() +
  geom_point(data = dat, aes(x = sqrt(prec_niche_breadth), y = log(area)), size = 2) + 
  geom_smooth(data = dat, aes(x = sqrt(prec_niche_breadth), y = log(area)), method = "lm", 
              se = FALSE, linewidth = 1, linetype = 1, colour = "black") +
  geom_smooth(data = dat[dat$Domiciliary == 1, ], aes(x = sqrt(prec_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 2, colour = "steelblue3") +
  geom_smooth(data = dat[dat$Peridomiciliary == 1, ], aes(x = sqrt(prec_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 1, colour = "orange2") +
  geom_smooth(data = dat[dat$Sylvatic == 1, ], aes(x = sqrt(prec_niche_breadth), y = log(area)), 
              method = "lm", se = FALSE, linewidth = 1, linetype = 1, colour = "chartreuse3") +
  labs(x = "Precipitation niche breadth", y = expression(paste("log Range size (km"^"2", ")"))) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none") +
  guides(colour = FALSE)

pdf("figures/Figure3.pdf", width = 22, height = 7)
wrap_plots(p1, p2, p3, ncol = 3)
dev.off()

# PGLS

comp_dat <- list()
for (i in 1:length(trees)) {
  comp_dat[[i]] <- comparative.data(data = dat, phy = trees[[i]],
                                    names.col = "species", vcv.dim = 2, 
                                    warn.dropped = TRUE)
}

pgls_spp <- pgls_gen <- pgls_fam <- pgls_ord <- pgls_cla <- list()
for (i in 1:length(comp_dat)) {
  pgls_spp[[i]] <- pgls(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) +
                        sqrt(prec_niche_breadth), lambda = "ML", data = comp_dat[[i]])
  pgls_gen[[i]] <- pgls(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) +
                        sqrt(prec_niche_breadth), lambda = "ML", data = comp_dat[[i]])
  pgls_fam[[i]] <- pgls(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) +
                        sqrt(prec_niche_breadth), lambda = "ML", data = comp_dat[[i]])
  pgls_ord[[i]] <- pgls(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) +
                        sqrt(prec_niche_breadth), lambda = "ML", data = comp_dat[[i]])
  pgls_cla[[i]] <- pgls(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) +
                        sqrt(prec_niche_breadth), lambda = "ML", data = comp_dat[[i]])
}

stats_pgls <- function(x) {
  
  df <- as.data.frame(matrix(ncol = 15, nrow = length(x)))
  colnames(df) <- c("intercept", "estimate_DBR", "Std_Error_DBR", "t_DBR", "p_DBR",
                    "estimate_temp", "Std_Error_temp", "t_temp", "p_temp",
                    "estimate_prec", "Std_Error_prec", "t_prec", "p_prec",
                    "R_squared", "Lambda")
  for (i in 1:length(x)) {
    df[i, 1] <- summary(x[[i]])$coefficients[1, 1]
    df[i, 2] <- summary(x[[i]])$coefficients[2, 1]
    df[i, 3] <- summary(x[[i]])$coefficients[2, 2]
    df[i, 4] <- summary(x[[i]])$coefficients[2, 3]
    df[i, 5] <- summary(x[[i]])$coefficients[2, 4]
    df[i, 6] <- summary(x[[i]])$coefficients[3, 1]
    df[i, 7] <- summary(x[[i]])$coefficients[3, 2]
    df[i, 8] <- summary(x[[i]])$coefficients[3, 3]
    df[i, 9] <- summary(x[[i]])$coefficients[3, 4]
    df[i, 10] <- summary(x[[i]])$coefficients[4, 1]
    df[i, 11] <- summary(x[[i]])$coefficients[4, 2]
    df[i, 12] <- summary(x[[i]])$coefficients[4, 3]
    df[i, 13] <- summary(x[[i]])$coefficients[4, 4]
    df[i, 14] <- summary(x[[i]])$r.squared
    df[i, 15] <- summary(x[[i]])$param[2]
  }
  
  return(df)
  
}

pgls_spp <- stats_pgls(pgls_spp)
pgls_gen <- stats_pgls(pgls_gen)
pgls_fam <- stats_pgls(pgls_fam)
pgls_ord <- stats_pgls(pgls_ord)
pgls_cla <- stats_pgls(pgls_cla)

res_pgls <- rbind(
  paste0(apply(pgls_spp, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_spp, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_spp, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_gen, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_gen, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_gen, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_fam, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_fam, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_fam, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_ord, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_ord, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_ord, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_cla, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_cla, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_cla, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")")
)

res_pgls <- as.data.frame(cbind(c("species", "genus", "family", "order", "class"), 
                                res_pgls))
colnames(res_pgls) <- c("DBR", colnames(pgls_spp))

write.csv(res_pgls, "tables/Table_pgls.csv")

## habitats

comp_dat_dom <- comp_dat_per <- comp_dat_syl <- list()
for (i in 1:length(trees)) {
  comp_dat_dom[[i]] <- comparative.data(data = dat[dat$Domiciliary == 1, ], 
                                        phy = keep.tip(trees[[i]], dat$species[dat$Domiciliary == 1]),
                                        names.col = "species", vcv.dim = 2, 
                                        warn.dropped = TRUE)
  comp_dat_per[[i]] <- comparative.data(data = dat[dat$Peridomiciliary == 1, ], 
                                        phy = keep.tip(trees[[i]], dat$species[dat$Peridomiciliary == 1]),
                                        names.col = "species", vcv.dim = 2, 
                                        warn.dropped = TRUE)
  comp_dat_syl[[i]] <- comparative.data(data = dat[dat$Sylvatic == 1, ], 
                                        phy = keep.tip(trees[[i]], dat$species[dat$Sylvatic == 1]),
                                        names.col = "species", vcv.dim = 2, 
                                        warn.dropped = TRUE)
}

pgls_spp_dom <- pgls_gen_dom <- pgls_fam_dom <- pgls_ord_dom <- pgls_cla_dom <- list()
for (i in 1:length(comp_dat_dom)) {
  pgls_spp_dom[[i]] <- pgls(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) +
                            sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_dom[[i]])
  pgls_gen_dom[[i]] <- pgls(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) +
                            sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_dom[[i]])
  pgls_fam_dom[[i]] <- pgls(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) +
                            sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_dom[[i]])
  pgls_ord_dom[[i]] <- pgls(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) +
                            sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_dom[[i]])
  pgls_cla_dom[[i]] <- pgls(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) +
                            sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_dom[[i]])
}

pgls_spp_per <- pgls_gen_per <- pgls_fam_per <- pgls_ord_per <- pgls_cla_per <- list()
for (i in 1:length(comp_dat_per)) {
  pgls_spp_per[[i]] <- pgls(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_per[[i]])
  pgls_gen_per[[i]] <- pgls(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_per[[i]])
  pgls_fam_per[[i]] <- pgls(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_per[[i]])
  pgls_ord_per[[i]] <- pgls(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_per[[i]])
  pgls_cla_per[[i]] <- pgls(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_per[[i]])
}

pgls_spp_syl <- pgls_gen_syl <- pgls_fam_syl <- pgls_ord_syl <- pgls_cla_syl <- list()
for (i in 1:length(comp_dat_syl)) {
  pgls_spp_syl[[i]] <- pgls(log(area) ~ log(speciesDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_syl[[i]])
  pgls_gen_syl[[i]] <- pgls(log(area) ~ log(genusDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_syl[[i]])
  pgls_fam_syl[[i]] <- pgls(log(area) ~ log(familyDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_syl[[i]])
  pgls_ord_syl[[i]] <- pgls(log(area) ~ log(orderDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_syl[[i]])
  pgls_cla_syl[[i]] <- pgls(log(area) ~ log(classDBR) + sqrt(temp_niche_breadth) +
                              sqrt(prec_niche_breadth), lambda = "ML", 
                            data = comp_dat_syl[[i]])
}

pgls_spp_dom <- stats_pgls(pgls_spp_dom)
pgls_gen_dom <- stats_pgls(pgls_gen_dom)
pgls_fam_dom <- stats_pgls(pgls_fam_dom)
pgls_ord_dom <- stats_pgls(pgls_ord_dom)
pgls_cla_dom <- stats_pgls(pgls_cla_dom)

pgls_spp_per <- stats_pgls(pgls_spp_per)
pgls_gen_per <- stats_pgls(pgls_gen_per)
pgls_fam_per <- stats_pgls(pgls_fam_per)
pgls_ord_per <- stats_pgls(pgls_ord_per)
pgls_cla_per <- stats_pgls(pgls_cla_per)

pgls_spp_syl <- stats_pgls(pgls_spp_syl)
pgls_gen_syl <- stats_pgls(pgls_gen_syl)
pgls_fam_syl <- stats_pgls(pgls_fam_syl)
pgls_ord_syl <- stats_pgls(pgls_ord_syl)
pgls_cla_syl <- stats_pgls(pgls_cla_syl)

res_pgls_habitat <- rbind(
  paste0(apply(pgls_spp_dom, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_spp_dom, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_spp_dom, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_gen_dom, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_gen_dom, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_gen_dom, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_fam_dom, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_fam_dom, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_fam_dom, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_ord_dom, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_ord_dom, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_ord_dom, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_cla_dom, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_cla_dom, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_cla_dom, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_spp_per, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_spp_per, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_spp_per, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_gen_per, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_gen_per, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_gen_per, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_fam_per, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_fam_per, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_fam_per, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_ord_per, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_ord_per, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_ord_per, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_cla_per, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_cla_per, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_cla_per, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_spp_syl, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_spp_syl, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_spp_syl, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_gen_syl, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_gen_syl, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_gen_syl, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_fam_syl, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_fam_syl, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_fam_syl, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_ord_syl, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_ord_syl, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_ord_syl, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")"),
  
  paste0(apply(pgls_cla_syl, 2, function(x) round(median(x), 3)), " (", 
         apply(pgls_cla_syl, 2, function(x) round(quantile(x, probs = 0.025), 3)), "-", 
         apply(pgls_cla_syl, 2, function(x) round(quantile(x, probs = 0.975), 3)), ")")
)

res_pgls_habitat <- as.data.frame(
  cbind(rep(c("Domiciliary", "Peridomiciliary", "Sylvatic"), each = 5),
        rep(c("species", "genus", "family", "order", "class"), 3),
        res_pgls_habitat)
  )
colnames(res_pgls_habitat) <- c("Habitat", "DBR", colnames(pgls_spp_dom))

write.csv(res_pgls_habitat, "tables/Table_pgls_habitat.csv")
