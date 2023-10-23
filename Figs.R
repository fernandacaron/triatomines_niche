rm(list = ls())

setwd("~/Documents/lab/triatomines")

library(phytools)
library(ggtree)
library(deeptime)
library(geiger)
library(viridis)
library(correlation)
library(phylosignal)
library(phylobase)
library(caper)
library(MuMIn)
library(wesanderson)
library(ggplot2)
library(MetBrewer)

tree <- read.nexus("data/tr_cal_f.tre")
trees <- read.nexus("data/1k_random_tria.trees")

dat <- read.delim("data/all_habitats_analisis.habit.txt", header = T)
row.names(dat) <- dat$species

#eliminate this samples, they are not in the phylogeny
dat <- dat[-c(24, 27), ]

name.check(tree, dat)

tree <- drop.tip(tree, setdiff(tree$tip.label, dat$species))
trees.pruned <- list()
for (i in 1:1000) {
	trees.pruned[[i]] <- drop.tip(trees[[i]], setdiff(trees[[i]]$tip.label,
	                                                  dat$species))
}
name.check(tree, dat)

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

# Analyses

#Verify if there is collinearity amog variables (OK)
dat_corr <- dat[, c("area", "niche_breadth_temp", "niche_breadth_precip", 
                    "PC1", "PC2", "PC3")]

result <- correlation(dat_corr)
s <- summary(result)
plot(s)

# Phylogenetic Signal
dat_red <- dat[, c("species", "classDBR", "orderDBR", "familyDBR", 
                   "genusDBR", "speciesDBR", "area", "niche_breadth_temp", 
                   "niche_breadth_precip", "PC1", "PC2", "PC3")]
dat_red <- dat_red[, -1]

# MCC
p4d_tr <- phylo4d(tree, dat_red)
phylosig_tr <- phyloSignal(p4d = p4d_tr, method = c("Lambda", "K"), 
                           reps = 1000)

# 1000 replicates
K <- K_p <- L <- L_p <- matrix(nrow = 11, ncol = 1000)
rownames(K) <- rownames(K_p) <- rownames(L) <- rownames(L_p) <- 
	colnames(dat_red)
for (i in 1:1000) {
	p4d <- phylo4d(trees.pruned[[i]], dat_red)
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

res_phylosig <- round(res_phylosig, 3)

write.csv(res_phylosig, "tables/Table1_unformatted.csv")

# PGLS

layout(matrix(1:6, ncol = 3))
hist(dat$niche_breadth_temp)
hist(dat$niche_breadth_precip)
hist(dat$area)
hist(dat$PC1)
hist(dat$PC2)
hist(dat$PC3)

# Scaled the predictors
dat$scale_nbtemp <- scale(dat$niche_breadth_temp, center = TRUE)
dat$scale_nbprec <- scale(dat$niche_breadth_precip, center = TRUE)
dat$scale_area <- scale(dat$area, center = TRUE)
dat$scale_pc1 <- scale(dat$PC1, center = TRUE)
dat$scale_pc2 <- scale(dat$PC2, center = TRUE)
dat$scale_pc3 <- scale(dat$PC3, center = TRUE)

layout(matrix(1:6, ncol = 3))
hist(dat$scale_nbtemp)
hist(dat$scale_nbprec)
hist(dat$scale_area)
hist(dat$scale_pc1)
hist(dat$scale_pc2)
hist(dat$scale_pc3)

col <- rgb(0.654902, 0.654902, 0.654902, 0.5)

col1 <- rgb(0.4470588, 0.5803922, 0.8313725)
col2 <- rgb(0.9921569, 0.3921569, 0.4039216)

col1_alp <- rgb(, 0.5803922, 0.8313725, 0.5)
col2_alp <- rgb(0.9921569, 0.3921569, 0.4039216, 0.5)

plot_pgls <- function(phy, data, formula, col1, col2, add) {

	comp_dat <- comparative.data(data = data, phy = phy,                         names.col = "species", vcv.dim = 2,
	                             warn.dropped = TRUE)

	res_pgls <- pgls(formula, lambda = "ML", data = comp_dat)

	coef <- res_pgls$model$coef[2:7, 1]
	names(coef) <- rownames(res_pgls$model$coef)[2:7]
	
	uci <- matrix(nrow = 6, ncol = 2)
	rownames(uci) <- names(coef)
	
	for (i in 1:6) {
		uci[i, 1] <- coef[i] - 1.96 * summary(res_pgls)$coefficients[1+i, 2]
		uci[i, 2] <- coef[i] + 1.96 * summary(res_pgls)$coefficients[1+i, 2]
	}

	mod_sel <- dredge(res_pgls)

	mod_avg <- model.avg(mod_sel, subset = delta < 2, fit = TRUE)   

	coefplot(coef, uci, shift = 0.3, pch = 16, col = col1, add = add,
	         labels = c("Range area", "PC1", "PC2", "PC3", "PBR", "TBR"),
	         xlim = c(-6, 16))
	plot(mod_avg, intercept = F, add = T, col = col2, pch = 17) 

}

# speciesDBR

pdf("figures/Fig2.pdf")

ggplot() + 
  geom_point(aes(x = dat$scale_area, y = dat$speciesDBR), 
                 color = "#C969A1") +
  geom_smooth(aes(x = dat$scale_area, y = dat$speciesDBR), method = lm,
              color = "#C969A1", fill = "#C969A1", alpha = 0.3) + 
  geom_point(aes(x = dat$scale_area, y = dat$genusDBR), 
                 color = "#EE8577") +
  geom_smooth(aes(x = dat$scale_area, y = dat$genusDBR), method = lm, 
                 color = "#EE8577", fill = "#EE8577", alpha = 0.3) + 
  geom_point(aes(x = dat$scale_area, y = dat$familyDBR),
             color = "#98AB76") +
  geom_smooth(aes(x = dat$scale_area, y = dat$familyDBR), method = lm,
              color = "#98AB76", fill = "#98AB76", alpha = 0.3) + 
  geom_point(aes(x = dat$scale_area, y = dat$orderDBR),
             color = "#C5DAF6") +
  geom_smooth(aes(x = dat$scale_area, y = dat$orderDBR), method = lm,
              color = "#C5DAF6", fill = "#C5DAF6", alpha = 0.3) + 
  geom_point(aes(x = dat$scale_area, y = dat$classDBR),
             color = "#DEC5DA") +
  geom_smooth(aes(x = dat$scale_area, y = dat$classDBR), method = lm,
              color = "#DEC5DA", fill = "#DEC5DA", alpha = 0.3) +
  xlab("Range area") +
  ylab("Diet breadth") + 
  theme_classic()

dev.off()

pdf("figures/Fig2_species_2.pdf")
plot_pgls(phy = trees.pruned[[1]], data = dat, col1 = col, col2 = col, add = F,
          formula = speciesDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
lapply(trees.pruned[2:1000], plot_pgls, data = dat, col1 = col, col2 = col, 
       add = T, formula = speciesDBR ~ scale_area + scale_pc1 + scale_pc2 + 
       scale_pc3 + scale_nbprec + scale_nbtemp)
plot_pgls(phy = tree, data = dat, col1 = col1, col2 = col2, add = T,
          formula = speciesDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
dev.off()

pdf("figures/Fig2_genus_2.pdf")
plot_pgls(phy = trees.pruned[[1]], data = dat, col1 = col, col2 = col, add = F,
          formula = genusDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
lapply(trees.pruned[2:1000], plot_pgls, data = dat, col1 = col, col2 = col, 
       add = T, formula = genusDBR ~ scale_area + scale_pc1 + scale_pc2 + 
       scale_pc3 + scale_nbprec + scale_nbtemp)
plot_pgls(phy = tree, data = dat, col1 = col1, col2 = col2, add = T,
          formula = genusDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
dev.off()

pdf("figures/Fig2_family_2.pdf")
plot_pgls(phy = trees.pruned[[1]], data = dat, col1 = col, col2 = col, add = F,
          formula = familyDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
lapply(trees.pruned[2:1000], plot_pgls, data = dat, col1 = col, col2 = col, 
       add = T, formula = familyDBR ~ scale_area + scale_pc1 + scale_pc2 + 
       scale_pc3 + scale_nbprec + scale_nbtemp)
plot_pgls(phy = tree, data = dat, col1 = col1, col2 = col2, add = T,
          formula = familyDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
dev.off()

pdf("figures/Fig2_order_2.pdf")
plot_pgls(phy = trees.pruned[[1]], data = dat, col1 = col, col2 = col, add = F,
          formula = orderDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
lapply(trees.pruned[2:1000], plot_pgls, data = dat, col1 = col, col2 = col, 
       add = T, formula = orderDBR ~ scale_area + scale_pc1 + scale_pc2 + 
       scale_pc3 + scale_nbprec + scale_nbtemp)
plot_pgls(phy = tree, data = dat, col1 = col1, col2 = col2, add = T,
          formula = orderDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
dev.off()

pdf("figures/Fig2_class_2.pdf")
plot_pgls(phy = trees.pruned[[1]], data = dat, col1 = col, col2 = col, add = F,
          formula = classDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
lapply(trees.pruned[2:1000], plot_pgls, data = dat, col1 = col, col2 = col, 
       add = T, formula = classDBR ~ scale_area + scale_pc1 + scale_pc2 + 
       scale_pc3 + scale_nbprec + scale_nbtemp)
plot_pgls(phy = tree, data = dat, col1 = col1, col2 = col2, add = T,
          formula = classDBR ~ scale_area + scale_pc1 + scale_pc2 + 
          scale_pc3 + scale_nbprec + scale_nbtemp)
dev.off()
