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
dat_corr <- dat[, c("speciesDBR", "niche_breadth_temp", "niche_breadth_precip", 
                    "PC1", "PC2", "PC3")]

result <- correlation(dat_corr)
s <- summary(result)
plot(s)


#Verify if there is collinearity amog variables (OK)
dat_corr <- dat[, c("genusDBR", "niche_breadth_temp", "niche_breadth_precip", 
                    "PC1", "PC2", "PC3")]

result <- correlation(dat_corr)
s <- summary(result)
plot(s)


#Verify if there is collinearity amog variables (OK)
dat_corr <- dat[, c("familyDBR", "niche_breadth_temp", "niche_breadth_precip", 
                    "PC1", "PC2", "PC3")]

result <- correlation(dat_corr)
s <- summary(result)
plot(s)


#Verify if there is collinearity amog variables (OK)
dat_corr <- dat[, c("orderDBR", "niche_breadth_temp", "niche_breadth_precip", 
                    "PC1", "PC2", "PC3")]

result <- correlation(dat_corr)
s <- summary(result)
plot(s)


#Verify if there is collinearity amog variables (OK)
dat_corr <- dat[, c("classDBR", "niche_breadth_temp", "niche_breadth_precip", 
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

######################
#FIGURE 2

run_pgls <- function(phy, data, formula) {

    comp_data <<- comparative.data(data = data, phy = phy,                         names.col = "species", vcv.dim = 2,
	                             warn.dropped = TRUE)
	res_pgls <- pgls(formula = formula, data = comp_data, lambda = "ML")

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

    res <- list()
    res$coef <- coef
    res$uci <- uci
    res$mod_avg <- mod_avg

    res
}

col <- rgb(0.654902, 0.654902, 0.654902, 0.5)

col1 <- rgb(0.4470588, 0.5803922, 0.8313725)
col2 <- rgb(0.9921569, 0.3921569, 0.4039216)

col1_alp <- rgb(0.5803922, 0.8313725, 0.5)
col2_alp <- rgb(0.9921569, 0.3921569, 0.4039216, 0.5)


#PGLS
comp_data <- NA
spp_pgls <- lapply(trees.pruned, run_pgls, data = dat, 
                   formula = scale_area ~ speciesDBR + scale_pc1 + scale_pc2 +
                   scale_pc3 + scale_nbprec + scale_nbtemp)
comp_data <- NA
spp_pgls_tr <- run_pgls(tree, data = dat,
                      formula = scale_area ~ speciesDBR + scale_pc1 + scale_pc2 + scale_pc3 + scale_nbprec + scale_nbtemp)

par(mfrow = c(2, 2))
comp_data_test <- comparative.data(data = dat, phy = tree, 
                                   names.col = "species", vcv.dim = 2,
                                   warn.dropped = TRUE)
test <- pgls(scale_area ~ speciesDBR + scale_pc1 + scale_pc2 + scale_pc3 + 
             scale_nbprec + scale_nbtemp, lambda = "ML", data = comp_data_test)
plot(test)

comp_data <- NA
gen_pgls <- lapply(trees.pruned, run_pgls, data = dat, 
                   formula = scale_area ~ genusDBR + scale_pc1 + scale_pc2 +
                   scale_pc3 + scale_nbprec + scale_nbtemp)
comp_data <- NA
gen_pgls_tr <- run_pgls(tree, data = dat,
                      formula = scale_area ~ genusDBR + scale_pc1 + scale_pc2 + scale_pc3 + scale_nbprec + scale_nbtemp)

par(mfrow = c(2, 2))
comp_data_test <- comparative.data(data = dat, phy = tree, 
                                   names.col = "species", vcv.dim = 2,
                                   warn.dropped = TRUE)
test <- pgls(scale_area ~ genusDBR + scale_pc1 + scale_pc2 + scale_pc3 + 
             scale_nbprec + scale_nbtemp, lambda = "ML", data = comp_data_test)
plot(test)

comp_data <- NA
fam_pgls <- lapply(trees.pruned, run_pgls, data = dat, 
                   formula = scale_area ~ familyDBR + scale_pc1 + scale_pc2 +
                   scale_pc3 + scale_nbprec + scale_nbtemp)
comp_data <- NA
fam_pgls_tr <- run_pgls(tree, data = dat,
                      formula = scale_area ~ familyDBR + scale_pc1 + scale_pc2 + scale_pc3 + scale_nbprec + scale_nbtemp)

par(mfrow = c(2, 2))
comp_data_test <- comparative.data(data = dat, phy = tree, 
                                   names.col = "species", vcv.dim = 2,
                                   warn.dropped = TRUE)
test <- pgls(scale_area ~ familyDBR + scale_pc1 + scale_pc2 + scale_pc3 + 
             scale_nbprec + scale_nbtemp, lambda = "ML", data = comp_data_test)
plot(test)

comp_data <- NA
ord_pgls <- lapply(trees.pruned, run_pgls, data = dat, 
                   formula = scale_area ~ orderDBR + scale_pc1 + scale_pc2 +
                   scale_pc3 + scale_nbprec + scale_nbtemp)
comp_data <- NA
ord_pgls_tr <- run_pgls(tree, data = dat,
                      formula = scale_area ~ orderDBR + scale_pc1 + scale_pc2 + scale_pc3 + scale_nbprec + scale_nbtemp)

par(mfrow = c(2, 2))
comp_data_test <- comparative.data(data = dat, phy = tree, 
                                   names.col = "species", vcv.dim = 2,
                                   warn.dropped = TRUE)
test <- pgls(scale_area ~ orderDBR + scale_pc1 + scale_pc2 + scale_pc3 + 
             scale_nbprec + scale_nbtemp, lambda = "ML", data = comp_data_test)
plot(test)

comp_data <- NA
cla_pgls <- lapply(trees.pruned, run_pgls, data = dat, 
                   formula = scale_area ~ classDBR + scale_pc1 + scale_pc2 +
                   scale_pc3 + scale_nbprec + scale_nbtemp)
comp_data <- NA
cla_pgls_tr <- run_pgls(tree, data = dat,
                      formula = scale_area ~ classDBR + scale_pc1 + scale_pc2 + scale_pc3 + scale_nbprec + scale_nbtemp)
par(mfrow = c(2, 2))
comp_data_test <- comparative.data(data = dat, phy = tree, 
                                   names.col = "species", vcv.dim = 2,
                                   warn.dropped = TRUE)
test <- pgls(scale_area ~ classDBR + scale_pc1 + scale_pc2 + scale_pc3 + 
             scale_nbprec + scale_nbtemp, lambda = "ML", data = comp_data_test)
plot(test)

col <- rgb(0.654902, 0.654902, 0.654902, 0.5)
cols1 <- c("#9F1B28", "#D95C1E", "#3C5B3A", "#4565AE", "#813B80")
cols2 <- c("#A64487", "#DC9A1B", "#718964", "#A9C1DB", "#CFB5D1")

# DBR
plot_pgls <- function(ntrees, pgls_obj, pgls_obj_tr, col, cols1, cols2,
                      add, shift1, shift2, parm) {
  coefplot(pgls_obj[[1]]$coef[parm], pgls_obj[[1]]$uci[parm, 1], 
           pgls_obj[[1]]$uci[parm, 2], shift = shift1, pch = 16, 
           col = col, xlim = c(-1.5, 1.5), ylim = c(-21, 30), add = add,
           labels = parm)
  try(
    coefplot(summary(pgls_obj[[1]]$mod_avg)$coefmat.full[parm, 1],
             confint(pgls_obj[[1]]$mod_avg, full = T, parm = parm)[1],
             confint(pgls_obj[[1]]$mod_avg, full = T, parm = parm)[2],
             shift = shift2, pch = 17, col = col, add = T)
  )
  
  for (i in 2:ntrees) {
    coefplot(pgls_obj[[i]]$coef[parm], pgls_obj[[i]]$uci[parm, 1], 
             pgls_obj[[i]]$uci[parm, 2], shift = shift1, pch = 16, 
             col = col, xlim = c(-6, 16), add = T)
    try(
      coefplot(summary(pgls_obj[[i]]$mod_avg)$coefmat.full[parm, 1],
               confint(pgls_obj[[i]]$mod_avg, full = T, parm = parm)[1],
               confint(pgls_obj[[i]]$mod_avg, full = T, parm = parm)[2],
               shift = shift2, pch = 17, col = col, add = T)
    )
  }
  
  coefplot(pgls_obj_tr$coef[parm], pgls_obj_tr$uci[parm, 1], 
           pgls_obj_tr$uci[parm, 2], shift = shift1, pch = 16, add = T,
           col = cols1, width = 0.3, dotcex = 1.5, lwd = 1.5, staplelwd = 1.5)
  try(
    coefplot(summary(pgls_obj_tr$mod_avg)$coefmat.full[parm, 1],
             confint(pgls_obj_tr$mod_avg, full = T, parm = parm)[1],
             confint(pgls_obj_tr$mod_avg, full = T, parm = parm)[2],
             shift = shift2, pch = 17, col = cols2, add = T, width = 0.3,
             dotcex = 1.5, lwd = 1.5, staplelwd = 1.5)
  )
}

pdf("figures/Fig3.pdf", height = 12, width = 14)

layout(matrix(1:2, ncol = 2))

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], F, 28, 27, "speciesDBR")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, 25, 24, "genusDBR")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, 22, 21, "familyDBR")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, 19, 18, "orderDBR")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, 16, 15, "classDBR")

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], T, 11, 10, "scale_nbtemp")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, 8, 7, "scale_nbtemp")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, 5, 4, "scale_nbtemp")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, 2, 1, "scale_nbtemp")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, -1, -2, "scale_nbtemp")

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], T, -6, -7, "scale_nbprec")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, -9, -10, "scale_nbprec")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, -12, -13, "scale_nbprec")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, -15, -16, "scale_nbprec")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, -18, -19, "scale_nbprec")

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], F, 28, 27, "scale_pc1")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, 25, 24, "scale_pc1")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, 22, 21, "scale_pc1")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, 19, 18, "scale_pc1")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, 16, 15, "scale_pc1")

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], T, 11, 10, "scale_pc2")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, 8, 7, "scale_pc2")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, 5, 4, "scale_pc2")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, 2, 1, "scale_pc2")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, -1, -2, "scale_pc2")

plot_pgls(3, spp_pgls, spp_pgls_tr, col, cols1[1], cols2[1], T, -6, -7, "scale_pc3")
plot_pgls(3, gen_pgls, gen_pgls_tr, col, cols1[2], cols2[2], T, -9, -10, "scale_pc3")
plot_pgls(3, fam_pgls, fam_pgls_tr, col, cols1[3], cols2[3], T, -12, -13, "scale_pc3")
plot_pgls(3, ord_pgls, ord_pgls_tr, col, cols1[4], cols2[4], T, -15, -16, "scale_pc3")
plot_pgls(3, cla_pgls, cla_pgls_tr, col, cols1[5], cols2[5], T, -18, -19, "scale_pc3")

dev.off()
