setwd("/Users/macbookpro/Dropbox/triatominae/triatominae_dbr/")
library(sp)
library(raster)
library(maps)
library(vegan)
library(ape)
library(geiger)
library(phytools)
library(subplex)

biocl <- raster::getData("worldclim",var="bio",res=2.5)
setwd("/Users/macbookpro/Dropbox/triatominae/wc2-5/")
bio1<-raster("bio1.bil")
bio2<-raster("bio2.bil")
bio3<-raster("bio3.bil")
bio4<-raster("bio4.bil")
bio5<-raster("bio5.bil")
bio6<-raster("bio6.bil")
bio7<-raster("bio7.bil")
bio8<-raster("bio8.bil")
bio9<-raster("bio9.bil")
bio10<-raster("bio10.bil")
bio11<-raster("bio11.bil")
bio12<-raster("bio12.bil")
bio13<-raster("bio13.bil")
bio14<-raster("bio14.bil")
bio15<-raster("bio15.bil")
bio16<-raster("bio16.bil")
bio17<-raster("bio17.bil")
bio18<-raster("bio18.bil")
bio19<-raster("bio19.bil")





#upload macro environmental variables
setwd("/Volumes/ELEMENTS/modalgem_mac/macro_enviro/wc2.0_30s_bio")
aet<-raster("aet2.tif")
ai<-raster("ai2.tif")
alpha<-raster("alpha2.tif")
etyr<-raster("etyr2.tif")

alpha.2.5 <- raster::resample(alpha, bio1, method = "bilinear")
aet.2.5 <- raster::resample(aet, bio1, method = "bilinear")
ai.2.5 <- raster::resample(ai, bio1, method = "bilinear")
etyr.2.5 <- raster::resample(etyr, bio1, method = "bilinear")

#npp
setwd("/Volumes/DANIELMAC2/npp/")
npp08.18 <- raster("npp08.18.tiff")


#stack all variables
clim.var <- stack(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,
                  bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19,
                  alpha.2.5,aet.2.5,ai.2.5,etyr.2.5,npp08.18)
                  

#################
#PCA
# PCs analyses were repeated for each type of habitat
setwd("/Users/macbookpro/Dropbox/triatominae/triatominae_dbr/")
dat<-read.table("triatof.all.habits.coords.txt", header = T) #
#dat<-unique(dat)

coords<-data.frame(x=dat$decimalLongitude, y=dat$decimalLatitude)
oc <- SpatialPoints(coords, proj4string = bio1@crs)

#First, let's have a look at the data
#map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
#points(coords, col="red", pch=16)

#It seems that the previous problematic points have been fixed.
clim<-extract(clim.var, oc)

res<-aggregate(x=clim, by=list(dat$species), FUN="mean", na.rm = TRUE)

rownames(res)<-res[,1]
res<-res[,-1]
PCA<-prcomp(res, center = TRUE, scale. = TRUE)

#the first three pcs were chosen
screeplot(PCA, bstick=TRUE)


summary(PCA, cutoff=0.001, loadings=TRUE) #it seems that we should retain the first three PCs
plot(PCA)


s <- summary(PCA)
s$importance

PCA$rotation
pcfactorloadings <- PCA$rotation[,1:3]

print(PCA$rotation[,1:3],digits = 1, cutoff=.3)
dim(PCA)
eixos.pca <- PCA$x

## Let's use the first five axis, which explain 95%
res$PC1 <- eixos.pca[,1]
res$PC2 <- eixos.pca[,2]
res$PC3 <- eixos.pca[,3]


#write.table(res, "triatopcafinal-domi.txt")


#############################
#CONTMAP AND DTT

#all habitats
pca.tri <- read.table("all_habitats_analisis.txt", header = TRUE)

row.names(pca.tri) <- pca.tri$species

pca.tri <- pca.tri[, c(
  "area", "niche_breadth_temp",
  "niche_breadth_precip", "PC1", "PC2", "PC3")]


#domiciliary
pca.tri.d <- read.delim("domiciliary_analisis.txt", header = TRUE)

row.names(pca.tri.d) <- pca.tri.d$species

pca.tri.d <- pca.tri.d[, c("area_km2", "niche_breadth_temp_bio5_bio6",
                            "niche_breadth_precip_bio16_bio17", "PC1", "PC2", "PC3")]

#domi-peri
tri.v.dp <- read.delim("Domiciliary_Peridomiciliary_analisis.txt", header = T)
row.names(tri.v.dp) <- tri.v.dp$species

tri.v_data.dp <- tri.v.dp[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                              "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]

#peridomiciliary
tri.v.p <- read.delim("peridomicilliary_ananlisis.txt", header = T)
row.names(tri.v.p) <- tri.v.p$species

tri.v.p <- tri.v.p[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                       #"niche_breadth_precip.bio16.bio17.", 
                       "PC1", "PC2", "PC3")]


#sylvatic
tri.v.s <- read.delim("sylvatic_analisys.txt", header = T)
row.names(tri.v.s) <- tri.v.s$species

tri.v_data.s <- tri.v.s[, c("area..km2.", "niche_breadth_temp.bio5.bio6.",
                            "niche_breadth_precip.bio16.bio17.", "PC1", "PC2", "PC3")]

#starting with the phylogenetic data
tr<-read.nexus("tr_cal_f.tre")


#tr$tip.label<-names[,2]

dat1<-treedata(tr, pca.tri)
dat1.d<-treedata(tr, pca.tri.d)
dat.dp <- treedata(tr, tri.v_data.dp)
dat.p <- treedata(tr, tri.v.p)
dat.s <- treedata(tr, tri.v_data.s)
#dat2<-treedata(tr, res)



#all habitats
contMap(dat1$phy, dat1$data[,2], cex=0.5)
contMap(dat1$phy, dat1$data[,3], cex=0.5)
contMap(dat1$phy, dat1$data[,4], cex=0.5)
contMap(dat1$phy, dat1$data[,5], cex=0.5)
contMap(dat1$phy, dat1$data[,6], cex=0.5)

dtt(dat1$phy, dat1$data[,2], nsim=1000, plot = TRUE)
dtt(dat1$phy, dat1$data[,3], nsim=1000, plot = TRUE)
dtt(dat1$phy, dat1$data[,4], nsim=1000, plot = TRUE)
dtt(dat1$phy, dat1$data[,5], nsim=1000, plot = TRUE)
dtt(dat1$phy, dat1$data[,6], nsim=1000, plot = TRUE)




#domi
contMap(dat1.d$phy, dat1.d$data[,2], cex=0.5)
contMap(dat1.d$phy, dat1.d$data[,3], cex=0.5)
contMap(dat1.d$phy, dat1.d$data[,4], cex=0.5)
contMap(dat1.d$phy, dat1.d$data[,5], cex=0.5)
contMap(dat1.d$phy, dat1.d$data[,6], cex=0.5)

dtt(dat1.d$phy, dat1.d$data[,2], nsim=1000, plot = TRUE)
dtt(dat1.d$phy, dat1.d$data[,3], nsim=1000, plot = TRUE)
dtt(dat1.d$phy, dat1.d$data[,4], nsim=1000, plot = TRUE)
dtt(dat1.d$phy, dat1.d$data[,5], nsim=1000, plot = TRUE)
dtt(dat1.d$phy, dat1.d$data[,6], nsim=1000, plot = TRUE)


#domi-peri
contMap(dat.dp$phy, dat.dp$data[,2], cex=0.5)
contMap(dat.dp$phy, dat.dp$data[,3], cex=0.5)
contMap(dat.dp$phy, dat.dp$data[,4], cex=0.5)
contMap(dat.dp$phy, dat.dp$data[,5], cex=0.5)
contMap(dat.dp$phy, dat.dp$data[,6], cex=0.5)

dtt(dat.dp$phy, dat.dp$data[,2], nsim=1000, plot = TRUE)
dtt(dat.dp$phy, dat.dp$data[,3], nsim=1000, plot = TRUE)
dtt(dat.dp$phy, dat.dp$data[,4], nsim=1000, plot = TRUE)
dtt(dat.dp$phy, dat.dp$data[,5], nsim=1000, plot = TRUE)
dtt(dat.dp$phy, dat.dp$data[,6], nsim=1000, plot = TRUE)

#peri

contMap(dat.p$phy, dat.p$data[,2], cex=0.5)
contMap(dat.p$phy, dat.p$data[,3], cex=0.5)
contMap(dat.p$phy, dat.p$data[,4], cex=0.5)
contMap(dat.p$phy, dat.p$data[,5], cex=0.5)


dtt(dat.p$phy, dat.p$data[,2], nsim=1000, plot = TRUE)
dtt(dat.p$phy, dat.p$data[,3], nsim=1000, plot = TRUE)
dtt(dat.p$phy, dat.p$data[,4], nsim=1000, plot = TRUE)
dtt(dat.p$phy, dat.p$data[,5], nsim=1000, plot = TRUE)

#sylvatic
contMap(dat.s$phy, dat.s$data[,2], cex=0.5)
contMap(dat.s$phy, dat.s$data[,4], cex=0.5)
contMap(dat.s$phy, dat.s$data[,5], cex=0.5)
contMap(dat.s$phy, dat.s$data[,6], cex=0.5)

dtt(dat.s$phy, dat.s$data[,2], nsim=1000, plot = TRUE)
dtt(dat.s$phy, dat.s$data[,3], nsim=1000, plot = TRUE)
dtt(dat.s$phy, dat.s$data[,4], nsim=1000, plot = TRUE)
dtt(dat.s$phy, dat.s$data[,5], nsim=1000, plot = TRUE)
