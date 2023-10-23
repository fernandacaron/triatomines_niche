library(ggplot2)
library(dplyr)
library(ggthemes)
library(GGally)
library(car)
library(ape)
library(phytools)
library(phangorn)
library(geiger)
library(picante)
library(caper)#PGLS
library(MuMIn)#model selection
library(ape)
library(paleotree)
library(geiger)
library(phylosignal)
library(adephylo)
library(phylobase)
library(nlme)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggraph)
library(ggiraphExtra)
library(pander)
library(knitr)
library(kableExtra)
library(treedata.table)
library(paleotree)
library(nlme)
library(caper)
library(phylosignal)
library(TreeTools)
library(AICcmodavg)
library(phylosignal)
library(adephylo)
library(jtools)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(mgcv)
library(tidymv)
library(mmodely)
library(phytools)
library(strap)

library(correlation)
library(see)
library(ggplot2)




setwd("/Users/macbookpro/Dropbox/triatominae/")

#tri.v <- read.delim("diet_breath_area_project_var.txt", header = T)
tri.v <- read.delim("all_habitats_analisis.txt", header = T)

#hist(tri.v$area)
#hist(tri.v$niche_breadth_temp)
#hist(tri.v$niche_breadth_precip)
#hist(tri.v$PC1)
#hist(tri.v$PC2)
#hist(tri.v$PC3)


#Load the all habitats data
tri.v_data <- tri.v[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                        "area", "niche_breadth_temp",
                              "niche_breadth_precip", "PC1", "PC2", "PC3")]


#Verify if there is collinearity amog variables
tri.v_data.corr <- tri.v_data[, c("area", "niche_breadth_temp",
                                  "niche_breadth_precip", "PC1", "PC2", "PC3")]

result <- correlation(tri.v_data.corr)
s <- summary(result)
plot(s)



#eliminate this samples, they are not in the phylogeny
tri.v_data <- tri.v_data[-c(24, 27), ]



#Phylogeny:
tree2 <- read.nexus("tr_cal_f.tre") # 


row.names(tri.v_data) <- tri.v_data$species
name.check(tree2, tri.v_data)
tree_tri <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v_data$species))
name.check(tree_tri, tri.v_data)
tri.v_data<-tri.v_data[,-1]

#####
#Phylogenetic Signal
#all habitats
p4d=phylo4d(tree_tri, tri.v_data)
phyloSignal(p4d = p4d, method = "all", reps = 10000)

#The results are in the Table 2 in the main text

####

#domiciliary
tri.v.d <- read.delim("domiciliary_analisis.txt", header = T)

tri.v_data.d <- tri.v.d[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                            "area_km2", "niche_breadth_temp_bio5_bio6",
                            "niche_breadth_precip_bio16_bio17", "PC1", "PC2", "PC3")]
tri.v_data.d.cor <- tri.v.d[, c(
                            "area_km2", "niche_breadth_temp_bio5_bio6",
                            "niche_breadth_precip_bio16_bio17", "PC1", "PC2", "PC3")]


#Analyze the collinearity
result.d <- correlation(tri.v_data.d.cor)
s.d <- summary(result.d)
plot(s.d)


tri.v_data.d <- tri.v_data.d[-c(14), ]

tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v_data.d) <- tri.v_data.d$species
name.check(tree2, tri.v_data.d)
tree_tri.d <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v_data.d$species))
name.check(tree_tri.d, tri.v_data.d)

tri.v_data.d<-tri.v_data.d[,-1]


#Phylogenetic signal for domiciliary

sig.phylo.d=phylo4d(tree_tri.d, tri.v_data.d)
signal.p.d = phyloSignal(p4d = sig.phylo.d, method = "all", reps = 10000)
signal.p.d

#####


#peridomiciliary

tri.v.p <- read.delim("peridomicilliary_ananlisis.txt", header = T)

tri.v_data.p <- tri.v.p[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                            "area..km2.", "niche_breadth_temp.bio5.bio6.",
                            #"niche_breadth_precip.bio16.bio17.", 
                            "PC1", "PC2", "PC3")]

tri.v_data.p.corr <- tri.v.p[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                            "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]

#Analyze the collinearity
result.p <- correlation(tri.v_data.p.corr)
s.p <- summary(result.p)
plot(s.p) #There is correlation between PreNBr and PC1, I chose PC1 over PreNBr



#tri.v_data.p <- tri.v_data.p[-c(14), ]

tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v_data.p) <- tri.v_data.p$species
name.check(tree2, tri.v_data.p)
tree_tri.p <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v_data.p$species))
name.check(tree_tri.p, tri.v_data.p)

tri.v_data.p<-tri.v_data.p[,-1]

#Phylogenetic signal for peridomiciliary

sig.phylo.p=phylo4d(tree_tri.p, tri.v_data.p)

signal.p.p = phyloSignal(p4d = sig.phylo.p, method = "all", reps = 10000)

signal.p.p


####


#silvatic

tri.v.s <- read.delim("sylvatic_analisys.txt", header = T)

tri.v_data.s.corr <- tri.v.s[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                            "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]


tri.v_data.s <- tri.v.s[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                            "area..km2.", "niche_breadth_temp.bio5.bio6.",
                                 #"niche_breadth_precip.bio16.bio17.", 
                            "PC1", "PC2", "PC3")]

#Analyze the collinearity
result.s <- correlation(tri.v_data.s.corr)
s.s <- summary(result.s)
plot(s.s) #There is correlation between PreNBr and PC1, I chose PC1 over PreNBr



tri.v_data.s <- tri.v_data.s[-c(19), ]

tree2 <- read.nexus("tr_cal_f.tre") # 
row.names(tri.v_data.s) <- tri.v_data.s$species
name.check(tree2, tri.v_data.s)
tree_tri.s <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v_data.s$species))
name.check(tree_tri.s, tri.v_data.s)

tri.v_data.s<-tri.v_data.s[,-1]

#phylogenetic signal for sylvatic

sig.phylo.s=phylo4d(tree_tri.s, tri.v_data.s)

signal.p.s = phyloSignal(p4d = sig.phylo.s, method = "all", reps = 10000)

signal.p.s

#####

#domiciliary-peridomiciliary
tri.v.dp <- read.delim("Domiciliary_Peridomiciliary_analisis.txt", header = T)

tri.v_data.dp.corr <- tri.v.dp[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                              "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]

tri.v_data.dp <- tri.v.dp[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                              "area..km2.", "niche_breadth_temp.bio5.bio6.",
                                   "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]
#Analyze the collinearity
result.dp <- correlation(tri.v_data.dp.corr)
s.dp <- summary(result.dp)
plot(s.dp)

tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v_data.dp) <- tri.v_data.dp$species
name.check(tree2, tri.v_data.dp)
tree_tri.dp <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v_data.dp$species))
name.check(tree_tri.dp, tri.v_data.dp)

tri.v_data.dp<-tri.v_data.dp[,-1]

#Phylogenetic signal for domi-peri

sig.phylo.dp=phylo4d(tree_tri.dp, tri.v_data.dp)

signal.p.dp = phyloSignal(p4d = sig.phylo.dp, method = "all", reps = 10000)

signal.p.dp


##################
#PGLS

#########################
####ALL HABITATS

tri.v <- read.delim("all_habitats_analisis.txt", header = T)

tri.v <- tri.v[-c(24, 27), ]

#Phylogeny:
tree2 <- read.nexus("tr_cal_f.tre") # 

#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v) <- tri.v$species
name.check(tree2, tri.v)
tree_tri <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v$species))
name.check(tree_tri, tri.v)
#tri.v_data<-tri.v_data[,-1]


#Scaled the predictors
tri.v$scale_nbtemp <- scale(tri.v$niche_breadth_temp, center = TRUE)
tri.v$scale_nbprec <- scale(tri.v$niche_breadth_precip, center = TRUE)
tri.v$scale_area <- scale(tri.v$area, center = TRUE)
tri.v$scale_pc1 <- scale(tri.v$PC1, center = TRUE)
tri.v$scale_pc2 <- scale(tri.v$PC2, center = TRUE)
tri.v$scale_pc3 <- scale(tri.v$PC3, center = TRUE)


tri.data <- comparative.data(data= tri.v, phy = tree_tri, 
                             names.col="species",  vcv.dim=2, warn.dropped=TRUE)  


#class diet breadth
class <- pgls(classDBR ~  scale_area * scale_pc1 + scale_nbprec   +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data)


summary(class)

#MODEL SELECTION
class_sel1 <- dredge(class)

#MODEL AVERAGE TAKING ACCOUNT MODELS WITH < 2 AICc   
class_avg1 <- model.avg(class_sel1, subset = delta < 2, fit = TRUE)    
summary(class_avg1)

###PLOT, THIS GRAPH COULD BE A GOOD MANNER TO SHOW OUR RESULTS TO AVOID THE TABLES
plot(class_avg1, full = NA, intercept = FALSE)


#RESIDUALS:
#Given that the averaged models represent several models where the coefficients 
#of each variable are estimated based on the weight of each independent model, we 
#checked the assumptions of normality and homogeneity of variance in the 
#distribution of residuals of a PGLS model containing all variables presented in 
#the averaged model 
par(mfrow=c(2,2))
class.resi <- pgls(classDBR ~  scale_area + scale_pc1 + scale_area*scale_pc1   +
                scale_nbprec  + 
                scale_pc2,
              lambda = "ML", data = tri.data)

plot(class.resi)


#plot class
#plot AREA

class.plot <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = mean(tri.v$scale_nbprec),
  scale_area = seq(min(tri.v$scale_area), max(tri.v$scale_area),length=180),
  scale_pc1 = mean(tri.v$scale_pc1),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

class.plot$pred <- predict(class_avg1, class.plot, full = TRUE)
plot_class.area <- ggplot(class.plot, aes(x = scale_area, y= pred)) +
  geom_point(data=tri.v, aes(x = scale_area, y = classDBR), size = 2, alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v,aes(x = scale_area, y = classDBR), method = "lm", se = T, col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
scale_y_continuous("Class diet breadth", limits = c(0, 6.5))+
  scale_x_continuous("Area range")


###
#AREA-PC1
class.plot.area.pc1 <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = mean(tri.v$scale_nbprec),
  #scale_area = seq(min(tri.v$scale_area), max(tri.v$scale_area),length=180),
  scale_area = (rep(c(-0.7913071, 3.60299),length=180)),
  scale_pc1 = seq(min(tri.v$scale_pc1), max(tri.v$scale_pc1),length=180),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))


min(class.plot.area.pc1$scale_area)
#-0.7913071
max(area.sort$scale_area)
#3.60299
class.plot.area.pc1$pred <- predict(class_avg1, class.plot.area.pc1, type = "response")

#class.pre_pc1$pred <- predict(class_avg, newdata = class.pre_pc1, type = "scale_pc1") 
#class.plot.area.pc1$scale_area1 <- class.pre_pc1$scale_area
class.plot.area.pc1$scale_area[class.plot.area.pc1$scale_area == -0.7913071] <- "Smaller Area"
#class.plot.area.pc1$scale_area1[class.plot.area.pc1$scale_area == 1.405841] <- "Media Area"
class.plot.area.pc1$scale_area[class.plot.area.pc1$scale_area == 3.60299] <- "Larger Area"

class.plot.area.pc1 %>% 
  ggplot() +
  aes(x = scale_pc1, y = pred, group = scale_area1, color = scale_area1) +
  geom_point(color = "grey", alpha = .7) +
  geom_smooth(method = "lm")


plot_class.area.pc1 <- ggplot(class.plot.area.pc1, aes(x = scale_pc1, y= pred, group = scale_area, color = scale_area)) +
  geom_point(data=tri.v, aes(x = scale_pc1, y = classDBR, color = scale_area), size = 2,alpha=0.6,col= "#2c7fb8")+
  #geom_smooth(data=tri.v,aes(x = scale_pc1, y = classDBR, color = scale_area), method = "lm", se = T)+
  geom_line(aes(lty = scale_area, color = scale_area), size = 1) +
  scale_linetype_manual(values=c("solid", "dashed")) + #2ca25f
  scale_color_manual(values=c("black",'black')) +
  #scale_colour_gradient(low = "#E2E0E0", high="#000000") +
                      #geom_line(aes(lty = scale_area), size = 1)+
  #scale_linetype_discrete(name='Area range', labels=c('Small','Media', 'Larger'))+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none")+
  #geom_line(linewidth = 1)+
  #scale_colour_colorblind()+
  scale_y_continuous("Class diet breadth", limits = c(0, 6.5))+
  scale_x_continuous("PC1")



#Order

order <- pgls(orderDBR ~scale_area * scale_pc1 + scale_nbprec   +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data)

order_sel <- dredge(order)

order_avg <- model.avg(order_sel, subset = delta < 2, fit = TRUE)    
summary(order_avg)  


#Residuals

order.resi <- pgls(orderDBR ~scale_area * scale_pc1 + scale_nbprec   +
                     scale_pc1  + 
                scale_area + scale_pc3,
              lambda = "ML", data = tri.data)

plot(order.resi)

#precipitation
order.plot <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = seq(min(tri.v$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v$scale_area),
  scale_pc1 = mean(tri.v$scale_pc1),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

order.plot$pred <- predict(order_avg, order.plot, full = TRUE)

plot_order.preci <- ggplot(order.plot, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v, aes(x = scale_nbprec, y = orderDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v,aes(x = scale_nbprec, y = orderDBR), method = "lm", se = T, col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none")+
  #geom_line(linewidth = 1)+
  
  #theme_modern()+
  #scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(1, 23.5))+
  scale_x_continuous("Precipitation niche breadth")


#area-pc1
#AREA-PC1
order.plot.area.pc1 <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = mean(tri.v$scale_nbprec),
  #scale_area = seq(min(tri.v$scale_area), max(tri.v$scale_area),length=130),
  scale_area = (rep(c(-0.7913071, 3.60299),length=180)),
  scale_pc1 = seq(min(tri.v$scale_pc1), max(tri.v$scale_pc1),length=180),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

#min(area.sort$scale_area)
#-0.7913071
#max(area.sort$scale_area)
#3.60299
order.plot.area.pc1$pred <- predict(order_avg, order.plot.area.pc1, type = "response")

#class.pre_pc1$pred <- predict(class_avg, newdata = class.pre_pc1, type = "scale_pc1") 
#order.plot.area.pc1$scale_area1 <- order.plot.area.pc1$scale_area
order.plot.area.pc1$scale_area[order.plot.area.pc1$scale_area == -0.7913071] <- "Smaller Area"
#order.plot.area.pc1$scale_area1[order.plot.area.pc1$scale_area == 0] <- "Media Area"
order.plot.area.pc1$scale_area[order.plot.area.pc1$scale_area == 3.60299] <- "Larger Area"





plot_order.area.pc1 <- ggplot(order.plot.area.pc1, aes(x = scale_pc1, y= pred,group = scale_area, color = scale_area)) +
  geom_point(data=tri.v, aes(x = scale_pc1, y = orderDBR, color = scale_area), size = 2,alpha=0.6, col= "#2c7fb8")+
  #geom_smooth(data=tri.v,aes(x = scale_pc1, y = orderDBR, color = scale_area), method = "lm", se = F)+
  geom_line(aes(lty = scale_area, color = scale_area), size = 1) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_color_manual(values=c("black",'black')) +  
#scale_colour_gradient(low = "#E2E0E0", high="#000000") +
  #geom_line(aes(lty = scale_area1), size = 1)+
  #scale_linetype_discrete(name='Area range', labels=c('Small','Media', 'Larger'))+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none"
  )+
  #geom_line(linewidth = 1)+
  #scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(1, 23.5))+
  scale_x_continuous("PC1")



#family
family <- pgls(familyDBR ~ scale_area * scale_pc1 + scale_nbprec   +
                 scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
               lambda = "ML", data = tri.data)

family_sel <- dredge(family)
family_avg <- model.avg(family_sel, subset = delta < 2, fit = TRUE)    
summary(family_avg)


#residuals
family.resi <- pgls(familyDBR ~ scale_area + scale_nbprec + scale_pc1 + scale_area*scale_pc1
                    +
                 scale_nbtemp,
               lambda = "ML", data = tri.data)

plot(family.resi)


#precipitation
family.plot <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = seq(min(tri.v$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v$scale_area),
  scale_pc1 = mean(tri.v$scale_pc1),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

family.plot$pred <- predict(family_avg, family.plot, full = TRUE)

plot_family.preci <- ggplot(family.plot, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v, aes(x = scale_nbprec, y = familyDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v,aes(x = scale_nbprec, y = familyDBR), method = "lm", se = T, col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Family diet breadth", limits = c(0, 45.5))+
  scale_x_continuous("Precipitation niche breadth")




#genus

genus <- pgls(genusDBR ~ scale_area * scale_pc1 + scale_nbprec   +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data)

#genus <- pgls(genusDBR ~ scale_area * scale_pc1 + scale_nbprec   *
# scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
# lambda = "ML", data = tri.data)

genus_sel <- dredge(genus)
summary(genus)
genus_avg <- model.avg(genus_sel, subset = delta < 2, fit = TRUE)    
summary(genus_avg)

#residuals
genus.resi <- pgls(genusDBR ~ scale_area * scale_pc1 + scale_area + scale_nbprec  +
                     scale_pc1 + scale_nbtemp,
              lambda = "ML", data = tri.data)

plot(genus.resi)

#precipitation
genus.plot <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = seq(min(tri.v$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v$scale_area),
  scale_pc1 = mean(tri.v$scale_pc1),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

genus.plot$pred <- predict(genus_avg, genus.plot, full = TRUE)

plot_genus.preci <- ggplot(genus.plot, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v, aes(x = scale_nbprec, y = genusDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v,aes(x = scale_nbprec, y = genusDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Genus diet breadth", limits = c(0, 56.5))+
  scale_x_continuous("Precipitation niche breadth")



#species
species <- pgls(speciesDBR ~ scale_area * scale_pc1 + scale_nbprec   +
                  scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data)

species_sel <- dredge(species)
#summary(species)
species_avg <- model.avg(species_sel, subset = delta < 2, fit = TRUE)    
summary(species_avg)    

#residuals
species.resi <- pgls(speciesDBR ~ scale_area * scale_pc1 + scale_area+ scale_nbprec +
                  scale_pc1,
                lambda = "ML", data = tri.data)


plot(species.resi)

#precipitation
species.plot <- data.frame(
  scale_nbtemp = mean(tri.v$scale_nbtemp),
  scale_nbprec = seq(min(tri.v$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v$scale_area),
  scale_pc1 = mean(tri.v$scale_pc1),
  scale_pc2 = mean(tri.v$scale_pc2),
  scale_pc3 = mean(tri.v$scale_pc3))

species.plot$pred <- predict(species_avg, species.plot, full = TRUE)

plot_species.preci <- ggplot(species.plot, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v, aes(x = scale_nbprec, y = genusDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v,aes(x = scale_nbprec, y = genusDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(0, 59.5))+
  scale_x_continuous("Precipitation niche breadth")


######################################




###sylvatic habitat

tri.v.s <- read.delim("sylvatic_analisys.txt", header = T)

tri.v.s = tri.v.s[-c(19),]

tri.v_data.s <- tri.v.s[, c("species" , "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                            "area..km2.", "niche_breadth_temp.bio5.bio6.",
                                 "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]


#Phylogeny:
tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v.s) <- tri.v.s$species
name.check(tree2, tri.v.s)
tree_tri.s <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v.s$species))
name.check(tree_tri.s, tri.v.s)


#niche breadth precip has corrleation with pc1
tri.v_data.s$scale_nbtemp <- scale(tri.v_data.s$niche_breadth_temp.bio5.bio6., center = TRUE)
tri.v_data.s$scale_nbprec <- scale(tri.v_data.s$niche_breadth_precip.bio16.bio17., center = TRUE)
tri.v_data.s$scale_area <- scale(tri.v_data.s$area..km2., center = TRUE)
tri.v_data.s$scale_pc1 <- scale(tri.v_data.s$PC1, center = TRUE)
tri.v_data.s$scale_pc2 <- scale(tri.v_data.s$PC2, center = TRUE)
tri.v_data.s$scale_pc3 <- scale(tri.v_data.s$PC3, center = TRUE)

tri.data.s <- comparative.data(data= tri.v_data.s, phy = tree_tri.s, 
                             names.col="species",  vcv.dim=2, warn.dropped=TRUE)  



#class
class.s <- pgls(classDBR ~  scale_area * scale_pc1    +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data.s)

class_sel1.s <- dredge(class.s)
summary(class.s)

class_avg1.s <- model.avg(class_sel1.s, subset = delta < 2, fit = TRUE)    
summary(class_avg1.s)


#residuals
class.s.resi <- pgls(classDBR ~  scale_area + scale_pc1    +scale_pc1 * scale_area +
                  scale_nbtemp,
                lambda = "ML", data = tri.data.s)

plot(class.s.resi)


#AREA-PC1
class.plot.area.pc1.s <- data.frame(
  scale_nbtemp = mean(tri.v_data.s$scale_nbtemp),
  scale_nbprec = mean(tri.v_data.s$scale_nbprec),
  #scale_area = seq(min(tri.v$scale_area), max(tri.v$scale_area),length=130),
  scale_area = (rep(c(-0.5405485,  4.590317),length=180)),
  scale_pc1 = seq(min(tri.v_data.s$scale_pc1), max(tri.v_data.s$scale_pc1),length=180),
  scale_pc2 = mean(tri.v_data.s$scale_pc2),
  scale_pc3 = mean(tri.v_data.s$scale_pc3))


min(tri.v_data.s$scale_area)
max(tri.v_data.s$scale_area)
median(tri.v_data.s$scale_area)
#area.sort= as.data.frame(tri.v$scale_area)
#area.sort=tri.v[order(tri.v$scale_area),]
#median(area.sort$scale_area)
#-0.4315619
#min(area.sort$scale_area)
#-0.7913071
#max(area.sort$scale_area)
#3.60299
class.plot.area.pc1.s$pred <- predict(class_avg1.s, class.plot.area.pc1.s, type = "response")

#class.pre_pc1$pred <- predict(class_avg, newdata = class.pre_pc1, type = "scale_pc1") 
#class.plot.area.pc1.s$scale_area1 <- class.plot.area.pc1.s$scale_area
class.plot.area.pc1.s$scale_area[class.plot.area.pc1.s$scale_area == -0.5405485] <- "Smaller Area"
#class.plot.area.pc1.s$scale_area[class.plot.area.pc1.s$scale_area == 0] <- "Media Area"
class.plot.area.pc1.s$scale_area[class.plot.area.pc1.s$scale_area == 4.590317] <- "Larger Area"



plot_class.area.pc1.s <- ggplot(class.plot.area.pc1.s, aes(x = scale_pc1, y= pred, group = scale_area, color = scale_area)) +
  geom_point(data=tri.v_data.s, aes(x = scale_pc1, y = classDBR, color = scale_area), size = 2,alpha=0.6,col= "#2c7fb8")+
  #geom_smooth(data=tri.v,aes(x = scale_pc1, y = orderDBR, color = scale_area), method = "lm", se = F)+
  #geom_line(aes(lty = scale_area), size = 1)
  geom_line(aes(lty = scale_area, color = scale_area), size = 1) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_color_manual(values=c("black",'black')) +  
  #scale_colour_gradient(low = "#E2E0E0", high="#000000") +
  #scale_linetype_discrete(name='Area range', labels=c('Small','Media', 'Larger'))+
  theme(#panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16)
        #legend.position = "none"
  )+
  #geom_line(linewidth = 1)+
  #scale_colour_colorblind()+
  scale_y_continuous("Class diet breadth", limits = c(1, 5.5))+
  scale_x_continuous("PC1")


#order

order.s <- pgls(orderDBR ~scale_area * scale_pc1    +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data.s)

order_sel.s <- dredge(order.s)
summary(order.s)
order_avg.s <- model.avg(order_sel.s, subset = delta < 2, fit = TRUE)    
summary(order_avg.s)  

#residuals
order.s.resi <- pgls(orderDBR ~scale_pc1 + scale_nbtemp + scale_pc3 + scale_area +
                       scale_area * scale_pc1    +
                       scale_pc2,
                lambda = "ML", data = tri.data.s)


plot(order.s.resi)

#family
family.s <- pgls(familyDBR ~ scale_area * scale_pc1  +
                 scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
               lambda = "ML", data = tri.data.s)

family_sel.s <- dredge(family.s)
summary(family.s)
family_avg.s <- model.avg(family_sel.s, subset = delta < 2, fit = TRUE)    
summary(family_avg.s)

#residuals
family.s.resi <- pgls(familyDBR ~ scale_nbtemp + scale_pc1 + 
                        scale_pc2 + scale_area + scale_pc3,
                 lambda = "ML", data = tri.data.s)


plot(family.s.resi)

#genus
genus.s <- pgls(genusDBR ~ scale_area * scale_pc1    +
                scale_nbtemp  + 
                scale_pc2 * scale_area + scale_area * scale_pc3,
              lambda = "ML", data = tri.data.s)

genus_sel.s <- dredge(genus.s)
summary(genus.s)
genus_avg.s <- model.avg(genus_sel.s, subset = delta < 2, fit = TRUE)    
summary(genus_avg.s)

#residuals
genus.s.resi <- pgls(genusDBR ~ scale_nbtemp + scale_pc1    +
                       scale_pc2  + scale_pc3,
                lambda = "ML", data = tri.data.s)

plot(genus.s.resi)

#species
#Note that there was correlation between precnbr and pc1
species.s <- pgls(speciesDBR ~ scale_area * scale_pc1    +
                  scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.s)


species_sel.s <- dredge(species.s)
summary(species.s)
species_avg.s <- model.avg(species_sel.s, subset = delta < 2, fit = TRUE)    
summary(species_avg.s)   


#resi
species.s.resi <- pgls(speciesDBR ~   scale_pc1    +
                    scale_nbtemp  + scale_pc2 + scale_pc3,
                  lambda = "ML", data = tri.data.s)

plot(species.s.resi)


#############################



####Domiciliary habitats

tri.v.d <- read.delim("domiciliary_analisis.txt", header = T)
tri.v.d = tri.v.d[-c(14),]
tri.v_data.d <- tri.v.d[, c("area_km2", "niche_breadth_temp_bio5_bio6",
                            "niche_breadth_precip_bio16_bio17", "PC1", "PC2", "PC3")]



#Phylogeny:
tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v.d) <- tri.v.d$species
name.check(tree2, tri.v.d)
tree_tri.d <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v.d$species))
name.check(tree_tri.d, tri.v.d)


#niche breadth precip has corrleation with pc1
tri.v.d$scale_nbtemp <- scale(tri.v.d$niche_breadth_temp_bio5_bio6, center = TRUE)
tri.v.d$scale_nbprec <- scale(tri.v.d$niche_breadth_precip_bio16_bio17, center = TRUE)
tri.v.d$scale_area <- scale(tri.v.d$area_km2, center = TRUE)
tri.v.d$scale_pc1 <- scale(tri.v.d$PC1, center = TRUE)
tri.v.d$scale_pc2 <- scale(tri.v.d$PC2, center = TRUE)
tri.v.d$scale_pc3 <- scale(tri.v.d$PC3, center = TRUE)

tri.data.d <- comparative.data(data= tri.v.d, phy = tree_tri.d, 
                               names.col="species",  vcv.dim=2, warn.dropped=TRUE)  



#class
class.d <- pgls(classDBR ~  scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.d)

class_sel1.d <- dredge(class.d)
#class_avg <- model.avg(class_sel, fit = TRUE)    
class_avg1.d <- model.avg(class_sel1.d, subset = delta < 2, fit = TRUE)    
summary(class_avg1.d)

#residuals
class.d.resi <- pgls(classDBR ~  scale_area * scale_pc1    + scale_nbprec +
                  scale_pc1  + scale_area,
                lambda = "ML", data = tri.data.d)


plot(class.d.resi)


#order

order.d <- pgls(orderDBR ~scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.d)

order_sel.d <- dredge(order.d)
order_avg.d <- model.avg(order_sel.d, subset = delta < 2, fit = TRUE)    
summary(order_avg.d)  

#residuals
order.d.resi <- pgls(orderDBR ~scale_area * scale_pc1 + 
                       scale_area   + scale_nbprec +
                  scale_pc1 + scale_pc2,
                lambda = "ML", data = tri.data.d)

plot(order.d.resi)

#precipitation
order.plot.d <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = seq(min(tri.v.d$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v.d$scale_area),
  scale_pc1 = mean(tri.v.d$scale_pc1),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))

order.plot.d$pred <- predict(order_avg.d, order.plot.d, full = TRUE)

plot_order.preci.d <- ggplot(order.plot.d, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v.d, aes(x = scale_nbprec, y = orderDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.d,aes(x = scale_nbprec, y = orderDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(0, 13.5))+
  scale_x_continuous("Precipitation niche breadth")


#family
family.d <- pgls(familyDBR ~ scale_area * scale_pc1    + scale_nbprec +
                   scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                 lambda = "ML", data = tri.data.d)

family_sel.d <- dredge(family.d)
family_avg.d <- model.avg(family_sel.d, subset = delta < 2, fit = TRUE)    
summary(family_avg.d)

#residuals:
family.d.resi <- pgls(familyDBR ~ scale_nbprec + scale_pc2 +scale_pc1+
                        scale_nbtemp + scale_area +
                        scale_area * scale_pc1,
                 lambda = "ML", data = tri.data.d)
  
plot(family.d.resi)

#precipitation
family.plot.d <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = seq(min(tri.v.d$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v.d$scale_area),
  scale_pc1 = mean(tri.v.d$scale_pc1),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))

family.plot.d$pred <- predict(family_avg.d, family.plot.d, full = TRUE)

plot_family.preci.d <- ggplot(family.plot.d, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v.d, aes(x = scale_nbprec, y = familyDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.d,aes(x = scale_nbprec, y = familyDBR), method = "lm", se = T, col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Family diet breadth", limits = c(0, 17.5))+
  scale_x_continuous("Precipitation niche breadth")


#genus
genus.d <- pgls(genusDBR ~ scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.d)

genus_sel.d <- dredge(genus.d)
genus_avg.d <- model.avg(genus_sel.d, subset = delta < 2, fit = TRUE)    
summary(genus_avg.d)

#resi
genus.d.resi <- pgls(genusDBR ~ scale_area + scale_nbprec +
                       scale_pc1  + scale_pc2 + 
                       scale_area *  scale_pc1 + scale_nbtemp + 
                       scale_area * scale_pc2,
                lambda = "ML", data = tri.data.d)

plot(genus.d.resi)

#precipitation
genus.plot.d <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = seq(min(tri.v.d$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v.d$scale_area),
  scale_pc1 = mean(tri.v.d$scale_pc1),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))

genus.plot.d$pred <- predict(genus_avg.d, genus.plot.d, full = TRUE)

plot_genus.preci.d <- ggplot(genus.plot.d, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v.d, aes(x = scale_nbprec, y = genusDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.d,aes(x = scale_nbprec, y = genusDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Genus diet breadth", limits = c(0, 20.5))+
  scale_x_continuous("Precipitation niche breadth")


#species
species.d <- pgls(speciesDBR ~ scale_area * scale_pc1    + scale_nbprec +
                    scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                  lambda = "ML", data = tri.data.d)

species_sel.d <- dredge(species.d)
species_avg.d <- model.avg(species_sel.d, subset = delta < 2, fit = TRUE)    
summary(species_avg.d) 

#resi
species.d.resi <- pgls(speciesDBR ~ scale_area + scale_nbprec +
                         scale_nbtemp +  scale_pc1    + scale_pc2 +
                    scale_pc1 * scale_area + scale_area * scale_pc2,
                  lambda = "ML", data = tri.data.d)


plot(species.d.resi)

#precipitation
species.plot.d <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = seq(min(tri.v.d$scale_nbprec), max(tri.v$scale_nbprec),length=180),
  scale_area = mean(tri.v.d$scale_area),
  scale_pc1 = mean(tri.v.d$scale_pc1),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))

species.plot.d$pred <- predict(species_avg.d, species.plot.d, full = TRUE)

plot_species.preci.d <- ggplot(species.plot.d, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v.d, aes(x = scale_nbprec, y = speciesDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.d,aes(x = scale_nbprec, y = speciesDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(0, 24.5))+
  scale_x_continuous("Precipitation niche breadth")

#PC1
species.plot.d.pc1 <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = mean(tri.v.d$scale_nbprec),
  scale_area = mean(tri.v.d$scale_area),
  scale_pc1 = seq(min(tri.v.d$scale_pc1), max(tri.v$scale_pc1),length=180),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))

species.plot.d.pc1$pred <- predict(species_avg.d, species.plot.d.pc1, full = TRUE)

plot_species.preci.d.pc1 <- ggplot(species.plot.d.pc1, aes(x = scale_pc1, y= pred)) +
  geom_point(data=tri.v.d, aes(x = scale_pc1, y = speciesDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.d,aes(x = scale_pc1, y = speciesDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  #scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(0, 24.5))+
  scale_x_continuous("PC1")


#AREA-PC1
species.plot.area.pc1.d <- data.frame(
  scale_nbtemp = mean(tri.v.d$scale_nbtemp),
  scale_nbprec = mean(tri.v.d$scale_nbprec),
  #scale_area = seq(min(tri.v$scale_area), max(tri.v$scale_area),length=130),
  scale_area = (rep(c(-0.593103, 3.05297),length=180)),
  scale_pc1 = seq(min(tri.v.d$scale_pc1), max(tri.v.d$scale_pc1),length=180),
  scale_pc2 = mean(tri.v.d$scale_pc2),
  scale_pc3 = mean(tri.v.d$scale_pc3))


min(tri.v.d$scale_area)
max(tri.v.d$scale_area)
median(tri.v.d$scale_area)
#area.sort= as.data.frame(tri.v$scale_area)
#area.sort=tri.v[order(tri.v$scale_area),]
#median(area.sort$scale_area)
#-0.4315619
#min(area.sort$scale_area)
#-0.7913071
#max(area.sort$scale_area)
#3.60299
species.plot.area.pc1.d$pred <- predict(species_avg.d, species.plot.area.pc1.d, type = "response")

#class.pre_pc1$pred <- predict(class_avg, newdata = class.pre_pc1, type = "scale_pc1") 
#species.plot.area.pc1.d$scale_area1 <- species.plot.area.pc1.d$scale_area
species.plot.area.pc1.d$scale_area[species.plot.area.pc1.d$scale_area == -0.593103] <- "Smaller Area"
#species.plot.area.pc1.d$scale_area1[species.plot.area.pc1.d$scale_area == 0] <- "Media Area"
species.plot.area.pc1.d$scale_area[species.plot.area.pc1.d$scale_area == 3.05297] <- "Larger Area"





plot_species.area.pc1.d <- ggplot(species.plot.area.pc1.d, aes(x = scale_pc1, y= pred, group = scale_area, color = scale_area)) +
  geom_point(data=tri.v.d, aes(x = scale_pc1, y = speciesDBR, color = scale_area), size = 2,alpha=0.6, col= "#2c7fb8")+
  #geom_smooth(data=tri.v,aes(x = scale_pc1, y = orderDBR, color = scale_area), method = "lm", se = F)+
  #geom_line(aes(lty = scale_area), size = 1)
  #scale_colour_gradient(low = "#E2E0E0", high="#000000") +
  geom_line(aes(lty = scale_area, color = scale_area), size = 1) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_color_manual(values=c("black",'black')) +  
    #scale_linetype_discrete(name='Area range', labels=c('Small','Media', 'Larger'))+
  theme(#panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.position = "none"
  )+
  #geom_line(linewidth = 1)+
  #scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(1, 24.5))+
  scale_x_continuous("PC1")

#########################


#Peridomiciliary habitats
#Phylogeny:
tree2 <- read.nexus("tr_cal_f.tre") # 


tri.v.p <- read.delim("peridomicilliary_ananlisis.txt", header = T)

tri.v.p <- tri.v.p[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                            "area..km2.", "niche_breadth_temp.bio5.bio6.",
                            #"niche_breadth_precip.bio16.bio17.", 
                            "PC1", "PC2", "PC3")]

#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v.p) <- tri.v.p$species
name.check(tree2, tri.v.p)
tree_tri.p <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v.p$species))
name.check(tree_tri.p, tri.v.p)


#niche breadth precip has corrleation with pc1
tri.v.p$scale_nbtemp <- scale(tri.v.p$niche_breadth_temp.bio5.bio6., center = TRUE)
#tri.v.p$scale_nbprec <- scale(tri.v.p$niche_breadth_precip.bio16.bio17., center = TRUE)
tri.v.p$scale_area <- scale(tri.v.p$area..km2., center = TRUE)
tri.v.p$scale_pc1 <- scale(tri.v.p$PC1, center = TRUE)
tri.v.p$scale_pc2 <- scale(tri.v.p$PC2, center = TRUE)
tri.v.p$scale_pc3 <- scale(tri.v.p$PC3, center = TRUE)

tri.data.p <- comparative.data(data= tri.v.p, phy = tree_tri.p, 
                               names.col="species",  vcv.dim=2, warn.dropped=TRUE)  



#class
#there was correlation between preecnbr and pc1
class.p <- pgls(classDBR ~  scale_area * scale_pc1 +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.p)

class_sel1.p <- dredge(class.p)
class_avg1.p <- model.avg(class_sel1.p, subset = delta < 2, 
                          fit = TRUE)    
summary(class_avg1.p)

#residuals
par(mfrow=c(2,2))
class.p.resi <- pgls(classDBR ~  scale_pc1 + scale_nbtemp +
                       scale_pc2 + scale_area + scale_area*scale_pc1 +
                  scale_pc2 * scale_area,
                lambda = "ML", data = tri.data.p)
plot(class.p.resi)

#order

order.p <- pgls(orderDBR ~scale_area * scale_pc1   +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.p)

order_sel.p <- dredge(order.p)
order_avg.p <- model.avg(order_sel.p, subset = delta < 2, fit = TRUE)    
summary(order_avg.p)  


#residuals
order.p.resi <- pgls(orderDBR ~ scale_pc1   +
                  scale_nbtemp  + 
                  scale_pc3,
                lambda = "ML", data = tri.data.p)

plot(order.p.resi)

class(tri.v.p)

#temp
order.plot.p.temp <- data.frame(
  scale_nbtemp = seq(min(tri.v.p$scale_nbtemp), max(tri.v.p$scale_nbtemp), length = 180),
  #scale_nbprec = mean(tri.v.p$scale_nbprec),
  scale_area = mean(tri.v.p$scale_area),
  scale_pc1 = mean(tri.v.p$scale_pc1),
  scale_pc2 = mean(tri.v.p$scale_pc2),
  scale_pc3 = mean(tri.v.p$scale_pc3))

order.plot.p.temp$pred <- predict(order_avg.p, order.plot.p.temp, full = TRUE)

plot_order.temp.p <- ggplot(order.plot.p.temp, aes(x = scale_nbtemp, y= pred)) +
  geom_point(data=tri.v.p, aes(x = scale_nbtemp, y = orderDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.p,aes(x = scale_nbtemp, y = orderDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(0, 17.5))+
  scale_x_continuous("Temperature niche breadth")

#pc1
order.plot.p.pc1 <- data.frame(
  scale_nbtemp = mean(tri.v.p$scale_nbtemp),
  #scale_nbprec = mean(tri.v.p$scale_nbprec),
  scale_area = mean(tri.v.p$scale_area),
  scale_pc1 = seq(min(tri.v.p$scale_pc1), max(tri.v.p$scale_pc1), length = 180),
  scale_pc2 = mean(tri.v.p$scale_pc2),
  scale_pc3 = mean(tri.v.p$scale_pc3))

order.plot.p.pc1$pred <- predict(order_avg.p, order.plot.p.pc1, full = TRUE)

plot_order.pc1.p <- ggplot(order.plot.p.pc1, aes(x = scale_pc1, y= pred)) +
  geom_point(data=tri.v.p, aes(x = scale_pc1, y = orderDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.p,aes(x = scale_pc1, y = orderDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(0, 17.5))+
  scale_x_continuous("PC1")

#family
family.p <- pgls(familyDBR ~ scale_area * scale_pc1    +
                   scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                 lambda = "ML", data = tri.data.p)

family_sel.p <- dredge(family.p)
summary(family_sel.p)
#Just one model
family_avg.p <- model.avg(family_sel.p, #subset = delta < 2, 
                          fit = TRUE)    
summary(family_avg.p)


#resi:
family.p.resi <- pgls(familyDBR ~ scale_nbtemp + scale_pc1    + scale_pc3 +
                        scale_pc2 + scale_area + scale_area * scale_pc1 + scale_pc2 * scale_area + 
                        scale_area * scale_pc3,
                 lambda = "ML", data = tri.data.p)  
  

plot(family.p.resi)

#temp
family.plot.p.temp <- data.frame(
  scale_nbtemp = seq(min(tri.v.p$scale_nbtemp), max(tri.v.p$scale_nbtemp), length = 180),
  #scale_nbprec = mean(tri.v.p$scale_nbprec),
  scale_area = mean(tri.v.p$scale_area),
  scale_pc1 = mean(tri.v.p$scale_pc1),
  scale_pc2 = mean(tri.v.p$scale_pc2),
  scale_pc3 = mean(tri.v.p$scale_pc3))

family.plot.p.temp$pred <- predict(family_avg.p, family.plot.p.temp, full = TRUE)

plot_family.temp.p <- ggplot(family.plot.p.temp, aes(x = scale_nbtemp, y= pred)) +
  geom_point(data=tri.v.p, aes(x = scale_nbtemp, y = familyDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.p,aes(x = scale_nbtemp, y = familyDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Family diet breadth", limits = c(0, 21.5))+
  scale_x_continuous("Temperature niche breadth")


#genus
genus.p <- pgls(genusDBR ~ scale_area * scale_pc1     +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.p)

genus_sel.p <- dredge(genus.p)
genus_avg.p <- model.avg(genus_sel.p, #subset = AIC < 2,
                         fit = TRUE)    
summary(genus_avg.p)
#residuals
genus.p.resi <- pgls(genusDBR ~ scale_nbtemp + scale_pc1     +
                       scale_pc2 + scale_pc3 + scale_area   + 
                  scale_pc1 * scale_area + scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.p)


plot(genus.p.resi)

#temp
genus.plot.p.temp <- data.frame(
  scale_nbtemp = seq(min(tri.v.p$scale_nbtemp), max(tri.v.p$scale_nbtemp), length = 180),
  #scale_nbprec = mean(tri.v.p$scale_nbprec),
  scale_area = mean(tri.v.p$scale_area),
  scale_pc1 = mean(tri.v.p$scale_pc1),
  scale_pc2 = mean(tri.v.p$scale_pc2),
  scale_pc3 = mean(tri.v.p$scale_pc3))

genus.plot.p.temp$pred <- predict(genus_avg.p, genus.plot.p.temp, full = TRUE)

plot_genus.temp.p <- ggplot(genus.plot.p.temp, aes(x = scale_nbtemp, y= pred)) +
  geom_point(data=tri.v.p, aes(x = scale_nbtemp, y = genusDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.p,aes(x = scale_nbtemp, y = genusDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Genus diet breadth", limits = c(0, 25.5))+
  scale_x_continuous("Temperature niche breadth")

#species
species.p <- pgls(speciesDBR ~ scale_area * scale_pc1  +
                    scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                  lambda = "ML", data = tri.data.p)

species_sel.p <- dredge(species.p)
summary(species.p)
species_avg.p <- model.avg(species_sel.p, subset = delta < 2, 
                           fit = TRUE)    
summary(species_avg.p) 

#resi:
species.p.resi <- pgls(speciesDBR ~ scale_nbtemp + scale_pc1  +  scale_pc2 +
                         scale_pc3 + scale_area + scale_pc2 * scale_area + scale_pc1 * scale_area +
                         scale_area * scale_pc3,
                  lambda = "ML", data = tri.data.p)

plot(species.p.resi)

#temp
species.plot.p.temp <- data.frame(
  scale_nbtemp = seq(min(tri.v.p$scale_nbtemp), max(tri.v.p$scale_nbtemp), length = 180),
  #scale_nbprec = mean(tri.v.p$scale_nbprec),
  scale_area = mean(tri.v.p$scale_area),
  scale_pc1 = mean(tri.v.p$scale_pc1),
  scale_pc2 = mean(tri.v.p$scale_pc2),
  scale_pc3 = mean(tri.v.p$scale_pc3))

species.plot.p.temp$pred <- predict(species_avg.p, species.plot.p.temp, full = TRUE)

plot_species.temp.p <- ggplot(species.plot.p.temp, aes(x = scale_nbtemp, y= pred)) +
  geom_point(data=tri.v.p, aes(x = scale_nbtemp, y = speciesDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v.p,aes(x = scale_nbtemp, y = speciesDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(0, 28.5))+
  scale_x_continuous("Temperature niche breadth")

###########################
#domi-peri

tree2 <- read.nexus("tr_cal_f.tre") # 
#read.nexus(file="/Users/macbookpro/Dropbox/triatominae/1K/1k_random_tria.trees") -> trees.1k
row.names(tri.v.dp) <- tri.v.dp$species
name.check(tree2, tri.v.dp)
tree_tri.dp <- drop.tip(tree2, setdiff(tree2$tip.label, tri.v.dp$species))
name.check(tree_tri.dp, tri.v.dp)

tri.v.dp <- read.delim("Domiciliary_Peridomiciliary_analisis.txt", header = T)

tri.v_data.dp <- tri.v.dp[, c("species", "classDBR", "orderDBR", "familyDBR", "genusDBR", "speciesDBR",
                              "area..km2.", "niche_breadth_temp.bio5.bio6.",
                            "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]


#niche breadth precip has corrleation with pc1
tri.v_data.dp$scale_nbtemp <- scale(tri.v_data.dp$niche_breadth_temp.bio5.bio6., center = TRUE)
tri.v_data.dp$scale_nbprec <- scale(tri.v_data.dp$niche_breadth_precip.bio16.bio17., center = TRUE)
tri.v_data.dp$scale_area <- scale(tri.v_data.dp$area..km2., center = TRUE)
tri.v_data.dp$scale_pc1 <- scale(tri.v_data.dp$PC1, center = TRUE)
tri.v_data.dp$scale_pc2 <- scale(tri.v_data.dp$PC2, center = TRUE)
tri.v_data.dp$scale_pc3 <- scale(tri.v_data.dp$PC3, center = TRUE)

tri.data.dp <- comparative.data(data= tri.v_data.dp, phy = tree_tri.dp, 
                               names.col="species",  vcv.dim=2, warn.dropped=TRUE)  



#class
class.dp <- pgls(classDBR ~  scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.dp)
summary(class.dp)
class_sel1.dp <- dredge(class.dp)
#class_avg <- model.avg(class_sel, fit = TRUE)    
class_avg1.dp <- model.avg(class_sel1.dp, subset = delta < 2, 
                          fit = TRUE)    
summary(class_avg1.dp)

#resi
class.dp.resi <- pgls(classDBR ~  scale_nbprec + scale_pc2 + scale_pc3 +
                        scale_nbtemp + scale_pc1,
                 lambda = "ML", data = tri.data.dp)

plot(class.dp.resi)

#precipitation
class.plot.dp <- data.frame(
  scale_nbtemp = mean(tri.v_data.dp$scale_nbtemp),
  scale_nbprec = seq(min(tri.v_data.dp$scale_nbprec), max(tri.v_data.dp$scale_nbprec),length=180),
  scale_area = mean(tri.v_data.dp$scale_area),
  scale_pc1 = mean(tri.v_data.dp$scale_pc1),
  scale_pc2 = mean(tri.v_data.dp$scale_pc2),
  scale_pc3 = mean(tri.v_data.dp$scale_pc3))

class.plot.dp$pred <- predict(class_avg1.dp, class.plot.dp, full = TRUE)

par(mfrow = c(2,2))
par(mfrow=c(2,2),las=1,bty="l",oma=c(1,1,1,1))
par(mfrow = c(2, 2), mar = rep(1, 4))


plot_class.preci.dp <- ggplot(class.plot.dp, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v_data.dp, aes(x = scale_nbprec, y = classDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v_data.dp,aes(x = scale_nbprec, y = classDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Class diet breadth", limits = c(0, 4.5))+
  scale_x_continuous("Precipitation niche breadth")

#order

order.dp <- pgls(orderDBR ~scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.dp)

order_sel.dp <- dredge(order.dp)
order_avg.dp <- model.avg(order_sel.dp, subset = delta < 2, fit = TRUE)    
summary(order_avg.dp)  

#resi
order.dp.resi <- pgls(orderDBR ~scale_nbprec + scale_pc2 + scale_nbtemp +
                        scale_pc1 + scale_pc3 + scale_area,
                 lambda = "ML", data = tri.data.dp)


plot(order.dp.resi)

#precipitation
order.plot.dp <- data.frame(
  scale_nbtemp = mean(tri.v_data.dp$scale_nbtemp),
  scale_nbprec = seq(min(tri.v_data.dp$scale_nbprec), max(tri.v_data.dp$scale_nbprec),length=180),
  scale_area = mean(tri.v_data.dp$scale_area),
  scale_pc1 = mean(tri.v_data.dp$scale_pc1),
  scale_pc2 = mean(tri.v_data.dp$scale_pc2),
  scale_pc3 = mean(tri.v_data.dp$scale_pc3))

order.plot.dp$pred <- predict(order_avg.dp, order.plot.dp, full = TRUE)

plot_order.preci.dp <- ggplot(order.plot.dp, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v_data.dp, aes(x = scale_nbprec, y = orderDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v_data.dp,aes(x = scale_nbprec, y = orderDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Order diet breadth", limits = c(0, 13.5))+
  scale_x_continuous("Precipitation niche breadth")


#family
family.dp <- pgls(familyDBR ~ scale_area * scale_pc1    + scale_nbprec +
                   scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                 lambda = "ML", data = tri.data.dp)

family_sel.dp <- dredge(family.dp)
summary(family.dp)
family_avg.dp <- model.avg(family_sel.dp, subset = delta < 2, fit = TRUE)    
summary(family_avg.dp)


#resi
family.dp.resi <- pgls(familyDBR ~ scale_nbprec + scale_nbtemp + scale_pc1 +
                         scale_pc2 + scale_area + scale_pc3,
                  lambda = "ML", data = tri.data.dp)

plot(family.dp.resi)


#precipitation
family.plot.dp <- data.frame(
  scale_nbtemp = mean(tri.v_data.dp$scale_nbtemp),
  scale_nbprec = seq(min(tri.v_data.dp$scale_nbprec), max(tri.v_data.dp$scale_nbprec),length=180),
  scale_area = mean(tri.v_data.dp$scale_area),
  scale_pc1 = mean(tri.v_data.dp$scale_pc1),
  scale_pc2 = mean(tri.v_data.dp$scale_pc2),
  scale_pc3 = mean(tri.v_data.dp$scale_pc3))

family.plot.dp$pred <- predict(family_avg.dp, family.plot.dp, full = TRUE)

plot_family.preci.dp <- ggplot(family.plot.dp, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v_data.dp, aes(x = scale_nbprec, y = familyDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v_data.dp,aes(x = scale_nbprec, y = familyDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Family diet breadth", limits = c(0, 19.5))+
  scale_x_continuous("Precipitation niche breadth")

#genus
genus.dp <- pgls(genusDBR ~ scale_area * scale_pc1    + scale_nbprec +
                  scale_nbtemp  + 
                  scale_pc2 * scale_area + scale_area * scale_pc3,
                lambda = "ML", data = tri.data.dp)

genus_sel.dp <- dredge(genus.dp)
genus_avg.dp <- model.avg(genus_sel.dp, subset = delta < 2, fit = TRUE)    
summary(genus_avg.dp)

#resi
genus.dp.resi <- pgls(genusDBR ~ scale_nbprec + scale_nbtemp +scale_pc1    + 
                        scale_pc2 +  scale_pc3,
                 lambda = "ML", data = tri.data.dp)
plot(genus.dp.resi)

#precipitation
genus.plot.dp <- data.frame(
  scale_nbtemp = mean(tri.v_data.dp$scale_nbtemp),
  scale_nbprec = seq(min(tri.v_data.dp$scale_nbprec), max(tri.v_data.dp$scale_nbprec),length=180),
  scale_area = mean(tri.v_data.dp$scale_area),
  scale_pc1 = mean(tri.v_data.dp$scale_pc1),
  scale_pc2 = mean(tri.v_data.dp$scale_pc2),
  scale_pc3 = mean(tri.v_data.dp$scale_pc3))

genus.plot.dp$pred <- predict(genus_avg.dp, genus.plot.dp, full = TRUE)

plot_genus.preci.dp <- ggplot(genus.plot.dp, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v_data.dp, aes(x = scale_nbprec, y = genusDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v_data.dp,aes(x = scale_nbprec, y = genusDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Genus diet breadth", limits = c(0, 22.5))+
  scale_x_continuous("Precipitation niche breadth")


#species
species.dp <- pgls(speciesDBR ~ scale_area * scale_pc1    + scale_nbprec +
                    scale_nbtemp  + scale_pc2 * scale_area + scale_area * scale_pc3,
                  lambda = "ML", data = tri.data.dp)

species_sel.dp <- dredge(species.dp)
species_avg.dp <- model.avg(species_sel.dp, subset = delta < 2, fit = TRUE)    
summary(species_avg.dp)

#resi
species.dp.resi <- pgls(speciesDBR ~ scale_nbprec +scale_nbtemp +
                          scale_area + scale_pc1,
                   lambda = "ML", data = tri.data.dp)

plot(species.dp.resi)

#precipitation
species.plot.dp <- data.frame(
  scale_nbtemp = mean(tri.v_data.dp$scale_nbtemp),
  scale_nbprec = seq(min(tri.v_data.dp$scale_nbprec), max(tri.v_data.dp$scale_nbprec),length=180),
  scale_area = mean(tri.v_data.dp$scale_area),
  scale_pc1 = mean(tri.v_data.dp$scale_pc1),
  scale_pc2 = mean(tri.v_data.dp$scale_pc2),
  scale_pc3 = mean(tri.v_data.dp$scale_pc3))

species.plot.dp$pred <- predict(species_avg.dp, species.plot.dp, full = TRUE)

plot_species.preci.dp <- ggplot(species.plot.dp, aes(x = scale_nbprec, y= pred)) +
  geom_point(data=tri.v_data.dp, aes(x = scale_nbprec, y = speciesDBR), size = 2,alpha=0.7, col= "#2c7fb8")+
  geom_smooth(data=tri.v_data.dp,aes(x = scale_nbprec, y = speciesDBR), method = "lm", se = T,col="#2ca25f")+
  theme(panel.background = element_rect(fill="white",colour="black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  #legend.position = "none")+
  #geom_line(linewidth = 1)+
  scale_colour_colorblind()+
  scale_y_continuous("Species diet breadth", limits = c(0, 22.5))+
  scale_x_continuous("Precipitation niche breadth")


###########################################################





#tri.data <- comparative.data(data= tri.v, phy = tree_tri, 
                            # names.col="species", vcv.dim=2, warn.dropped=TRUE)           


####analysis linear
#nbtemp
plot(tri.v$niche_breadth_temp, tri.v$classDBR)
abline(lm(tri.v$classDBR ~ tri.v$niche_breadth_temp))
cl.temp <- lm(tri.v$classDBR ~ tri.v$niche_breadth_temp)
summary(cl.temp)

#nbpre
plot(tri.v$niche_breadth_precip_bio16_bio17, tri.v$classDBR)
abline(lm(tri.v$classDBR ~ tri.v$niche_breadth_precip_bio16_bio17))
cl.pre <- lm(tri.v$classDBR ~ tri.v$niche_breadth_precip)
summary(cl.pre)

#pc1
plot(tri.v$PC1, tri.v$classDBR)
abline(lm(tri.v$classDBR ~ tri.v$PC1))
cl.pc1 <- lm(tri.v$classDBR ~ tri.v$PC1 * tri.v$area + tri.v$niche_breadth_precip_bio16_bio17 * tri.v$area)
summary(cl.pc1)

#pc2
plot(tri.v$PC2, tri.v$classDBR)
abline(lm(tri.v$classDBR ~ tri.v$PC2))
cl.pc2 <- lm(tri.v$classDBR ~ tri.v$PC2 * tri.v$area)
summary(cl.pc2)

#pc3
plot(tri.v$PC3, tri.v$classDBR)
abline(lm(tri.v$classDBR ~ tri.v$PC3))
cl.pc3 <- lm(tri.v$classDBR ~ tri.v$PC3 * tri.v$area)
summary(cl.pc3)

#order
plot(tri.v$area, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$area))
or.ar <- lm(tri.v$orderDBR ~ tri.v$area)
summary(or.ar)

plot(tri.v$niche_breadth_temp_bio5_bio6, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$niche_breadth_temp_bio5_bio6))
or.temp <- lm(tri.v$orderDBR ~ tri.v$niche_breadth_temp_bio5_bio6)
summary(or.temp)

plot(tri.v$niche_breadth_precip_bio16_bio17, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$niche_breadth_precip_bio16_bio17))
or.pre <- lm(tri.v$orderDBR ~ tri.v$niche_breadth_precip_bio16_bio17)
summary(or.pre)

plot(tri.v$PC1, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$PC1))
or.pc1 <- lm(tri.v$orderDBR ~ tri.v$PC1)
summary(or.pc1)

plot(tri.v$PC2, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$PC2))
or.pc2 <- lm(tri.v$orderDBR ~ tri.v$PC2)
summary(or.pc2)

plot(tri.v$PC3, tri.v$orderDBR)
abline(lm(tri.v$orderDBR ~ tri.v$PC3))
or.pc3 <- lm(tri.v$orderDBR ~ tri.v$PC3)
summary(or.pc3)

#FAMILY
plot(tri.v$area, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$area))
fa.ar <- lm(tri.v$familyDBR ~ tri.v$area)
summary(fa.ar)

plot(tri.v$niche_breadth_temp_bio5_bio6, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$niche_breadth_temp_bio5_bio6))
fa.temp <- lm(tri.v$familyDBR ~ tri.v$niche_breadth_temp_bio5_bio6)
summary(fa.temp)

plot(tri.v$niche_breadth_precip_bio16_bio17, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$niche_breadth_precip_bio16_bio17))
fa.pre <- lm(tri.v$familyDBR ~ tri.v$niche_breadth_precip_bio16_bio17)
summary(fa.pre)

plot(tri.v$PC1, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$PC1))
fa.pc1 <- lm(tri.v$familyDBR ~ tri.v$PC1)
summary(fa.pc1)

plot(tri.v$PC2, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$PC2))
fa.pc2 <- lm(tri.v$familyDBR ~ tri.v$PC2)
summary(fa.pc2)

plot(tri.v$PC3, tri.v$familyDBR)
abline(lm(tri.v$familyDBR ~ tri.v$PC3))
fa.pc3 <- lm(tri.v$familyDBR ~ tri.v$PC3)
summary(fa.pc3)

#genus

plot(tri.v$area, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$area))
ge.ar <- lm(tri.v$genusDBR ~ tri.v$area)
summary(ge.ar)

plot(tri.v$niche_breadth_temp_bio5_bio6, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$niche_breadth_temp_bio5_bio6))
ge.temp <- lm(tri.v$genusDBR ~ tri.v$niche_breadth_temp_bio5_bio6)
summary(ge.temp)

plot(tri.v$niche_breadth_precip_bio16_bio17, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$niche_breadth_precip_bio16_bio17))
ge.pre <- lm(tri.v$genusDBR ~ tri.v$niche_breadth_precip_bio16_bio17)
summary(ge.pre)

plot(tri.v$PC1, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$PC1))
ge.pc1 <- lm(tri.v$genusDBR ~ tri.v$PC1)
summary(ge.pc1)

plot(tri.v$PC2, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$PC2))
ge.pc2 <- lm(tri.v$genusDBR ~ tri.v$PC2)
summary(ge.pc2)

plot(tri.v$PC3, tri.v$genusDBR)
abline(lm(tri.v$genusDBR ~ tri.v$PC3))
ge.pc3 <- lm(tri.v$genusDBR ~ tri.v$PC3)
summary(ge.pc3)


#species

plot(tri.v$area, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$area))
s.ar <- lm(tri.v$speciesDBR ~ tri.v$area)
summary(s.ar)

plot(tri.v$niche_breadth_temp_bio5_bio6, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$niche_breadth_temp_bio5_bio6))
s.temp <- lm(tri.v$speciesDBR ~ tri.v$niche_breadth_temp_bio5_bio6)
summary(s.temp)

plot(tri.v$niche_breadth_precip_bio16_bio17, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$niche_breadth_precip_bio16_bio17))
s.pre <- lm(tri.v$speciesDBR ~ tri.v$niche_breadth_precip_bio16_bio17)
summary(s.pre)

plot(tri.v$PC1, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$PC1))
s.pc1 <- lm(tri.v$speciesDBR ~ tri.v$PC1)
summary(s.pc1)

plot(tri.v$PC2, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$PC2))
s.pc2 <- lm(tri.v$speciesDBR ~ tri.v$PC2)
summary(s.pc2)

plot(tri.v$PC3, tri.v$speciesDBR)
abline(lm(tri.v$speciesDBR ~ tri.v$PC3))
s.pc3 <- lm(tri.v$speciesDBR ~ tri.v$PC3)
summary(s.pc3)

##########

######
#geoscale
### Tree in geological scales
tree <- read.nexus("tr_cal_f.tre")
#Find root for plotting first by checking ages of nodes (nodeHeights)
lengths <- nodeHeights(tree)

#Now let's find the biggest number, that's our root node
root.time <- max(lengths)

#Set root for plotting
tree$root.time <- root.time

# grab our OTUs (Operational Taxonomic Units = taxa)
all_otus <- tree$tip.label

# Create an empty matrix containing the taxa, this is required by strap
all_otudates <- matrix(0, nrow = length(all_otus), ncol=2)

# Turn the matrix into a data frame
all_otudates <- data.frame(all_otudates)

#set the row names to the taxa (OTUs)
row.names(all_otudates) <- all_otus

# set column names to FAD (First Appearance Datum) and LAD (Last Appearance Datum)
colnames(all_otudates) <- c('FAD','LAD')    

geoscalePhylo(tree,ages=all_otudates, cex.tip=0.5, cex.age=0.7, cex.ts=0.7, lwd=3, width=1, x.lim=c(-15,44), units=c("Period", "Epoch"), boxes="Epoch")


##################

