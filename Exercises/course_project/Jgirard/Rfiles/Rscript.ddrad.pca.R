###Taking output from Bluewaves and attempting to R data visualizations on my personal PC in R

setwd("~/Documents/Documents/URI/Ph.D/Disertation_Classes_TA/Course Work/Spring 2022/Bio594_SpTopics_Puritz/594mdfiles/Final.Project/FinalProject_R")

install.packages ("adegenet") 
install.packages ("vcfR")
install.packages ("hierfstat")
# loading packages
library(adegenet)
library(vcfR)

#creating objects from data files 

my_vcf <- read.vcfR("./SNP.TRSdp5p05FHWEmaf05.2A.recode.vcf")


my_genind <- vcfR2genind(my_vcf) #making a gene ID file from original VCF
#str (my_genind)# returns ..@strata : NULL
#region ------------------------------------------------------------
#Ask about this block of code.

#"strata" needs to be created. Right now it appears to come from nothing and is filled with my_geneID

strata<- read.table("./ddrad.popmap.tab.txt", header=FALSE) #What is strata? It is not in the directory. "LibraryInfo"?
strata_df <- data.frame(strata)#turn to data frame


strata(my_genind) <- strata_df #set strata variable in my_genind to the strata dataframe
#error: number of idv is not equal between the data frame and the my_genid file. 
#one of the individuals was removed during the adegenet process this can be round in the first column of the @tab matrix in genid
 my_genind@tab[,1]
 idvremaining <- my_genind@tab[,1]
 idvremaining.dat <- as.data.frame(idvremaining)
 retainedidv <- row.names(idvremaining.dat)
 retained.dat <- as.data.frame(retainedidv)
 write.csv(retained.dat, "Retained.idv.csv")

# Used the match function in excel to find the missing individual – R18_R18
# R18_R18 is located in row 130
# Using a new popmap file with row 130 removed from the original ddrad.popmap.tab.txt

ddrad.popmap.sub <- strata_df[-130,]
strata_df <- ddrad.popmap.sub
strata_df <-`colnames<-`(strata_df, c("ID","Population"))# giving column names my_genind will recognize. "Population" is important

strata(my_genind) <- strata_df #set strata variable in my_genind to the strata dataframe
# worked!

setPop(my_genind) <- ~Population #worked
#------------------------------------------------------------------------------ 

#endregion-------------------------------------------------------------

#Test Population Structure
library(hierfstat)
#fstat(my_genind)

# heirfstat documentation https://cran.r-project.org/web/packages/hierfstat/hierfstat.pdf
basic.stats(my_genind) #overall Fst here is 0.0788...very low popn structure


#pp.fst(my_genind)#should report fst per pair | did not run
#matFst <- pairwise.fst(my_genind)

#PCA building

X <- tab(my_genind, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50)) #Creating barplot of eigenvectors (PCA axes?)
s.class(pca1$li, pop(my_genind))
title("PCA E.Flexuosa populations across the Caribbean Basin\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- funky(15)
s.class(pca1$li, pop(my_genind),xax=1,yax=2, col=col, axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)# plots PCA

#Let's try ggplot I have a little more experience here

install.packages("ggplot2")
install.packages("ggfortify")
library(ggplot2)
library(ggfortify)

autoplot(pca1)# saddly autoplot does not play nice with this object type
# I should be able to pull out the points manually
# documentation says pca1@li are the principle component scores

pca1$li

pca1.2.dat <- pca1$li[,1:3] #take cols 1 - 3
pca1.2.dat # lookgs good

#Now to add the population data
#Thankfully data should already be in the same order as the popmap

pca1.2.dat$Population <- strata_df$Population
pca1.2.dat# looks good
#bringing in a table with the depth lineage variable
ddrad.popmap.depth <- read.table("ddrad.popmap.depth.tab.txt", header = FALSE, sep = "",
                                 col.names = c("ID", "Population","Depth"))
# Removing sample R18_R18
ddrad.popmap.depth.sub <- ddrad.popmap.depth[-130,]

# Adding depth lineage col to pca1.2 data set

pca1.2.dat$Depth <- ddrad.popmap.depth.sub$Depth
pca1.2.dat


install.packages("viridis")# color pallets for R many robust to colorblindness
library(viridis)

#By Population
ggplot(pca1.2.dat, aes(x = Axis1, y = Axis2, col=Population)) + geom_point(size = 2)+
  scale_color_viridis_d()+
  theme_classic()+ labs(x="PC1", y="PC2")

#By Depth
ggplot(pca1.2.dat, aes(x = Axis1, y = Axis2, col=Depth)) + geom_point(size = 2)+
  #scale_color_manual(aesthetics = c("Green", "Red", "Blue"))+
  scale_color_viridis_d()+
  theme_classic()+ labs(x="PC1", y="PC2")


##Attempting an admixture plot source (https://rpubs.com/tlama/scatterpie)
# Install packages
# Many of these I may not need since I do not plan on attaching these data to maps, yet.
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap) 
#my machine needs BiocManager to install "LEA" package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LEA")
library(LEA)

library(maps)
library(mapdata)
library(maptools)
library(mapproj)
library(mapplots)
library(dismo)
library(plotrix)
library(dplyr)
library(rgdal)
library(ggplot2)
library(raster)
library(TESS)
library(ggforce)
install.packages('scatterpie')
library(scatterpie)

#Read in VCF
my_vcf <- read.vcfR("./SNP.TRSdp5p05FHWEmaf05.2A.recode.vcf")
head(my_vcf)
#Designate input file
input.file <- "./SNP.TRSdp5p05FHWEmaf05.2A.recode.vcf"

str(extract.gt(my_vcf, element = 'DP', as.numeric=TRUE))
dp<- extract.gt(my_vcf, element = 'DP', as.numeric=TRUE)##why was this done in the referene code. It never comes back and below the line it says they removed one of their samples


vcf2geno(input.file = "SNP.TRSdp5p05FHWEmaf05.2A.recode.vcf", output.file = "SNP.TRSdp5p05HWEFmaf05.2A.recode.vcf.geno", force = TRUE)
#Error in vcf2geno(input.file = "./SNP.TRSdp5p05FHWEmaf05.2A.recode.vcf",  : 
#internal error in trio library

###No clear documentation on how to fix this. After much trying I cannot figure out how to "look" at the vcf file in a way that helps me to solve the problem

#ce <- cross.entropy(my_vcf, K= 2)
#best <- which.min(ce)
#qmatrixK2 = Q(my_vcf, K = 2, run = best)
#barplot(t(qmatrixK2))
