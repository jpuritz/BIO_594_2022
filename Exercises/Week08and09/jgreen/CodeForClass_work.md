# Code for Week 8 & Week 9 assignment

library(adegenet)

library(vcfR)

my_vcf <- read.vcfR("SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf")

my_genind <- vcfR2genind(my_vcf)

strata<- read.table("LibraryInfo", header=TRUE)

strata_df <- data.frame(strata)

strata(my_genind) <- strata_df

setPop(my_genind) <- ~Population

#Test Population Structure

library("hierfstat")

## fstat(my_genind)

summary(my_genind)

names(my_genind)

pop(my_genind)

basic.stats(my_genind)

#PCA

X <- tab(my_genind, freq = TRUE, NA.method = "mean")

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

s.class(pca1$li, pop(my_genind))

title("PCA of real dataset\naxes 1-2")

add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- funky(15)

s.class(pca1$li, pop(my_genind),xax=1,yax=2, col=col, axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)

#DAPC

grp <- find.clusters(my_genind)

table(pop(my_genind), grp$grp)

table.value(table(pop(my_genind), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("ori", 1:16))

dapc1 <- dapc(my_genind, grp$grp)

scatter(dapc1,col=col,bg="white", solid=1)

contrib <- loadingplot(dapc1$var.contr, axis=1, thres=.01, lab.jitter=1)

contrib

setPop(my_genind) <- ~Library

dapc1 <- dapc(my_genind, pop(my_genind))

contrib <- loadingplot(dapc1$var.contr, axis=1, thres=.05, lab.jitter=1)

#Structure Like

compoplot(dapc1, posi="bottomright",txt.leg=paste("Cluster", 1:2), lab="", ncol=1, xlab="individuals")

temp <- which(apply(dapc1$posterior,2, function(e) all(e<0.9)))

compoplot(dapc1, subset=temp, posi="bottomright", txt.leg=paste("Cluster", 1:2), ncol=2)

