source("plot_R.r")
plot_bayescan("SNP.TRSdp5p05FH_fst.txt")

#Load pcadapt library
#install.packages("pcadapt", dependencies = TRUE)
library(pcadapt)

#load our VCF file into R
filename <- read.pcadapt("SNP.TRSdp5p05FHWE2A.recode.vcf", type = "vcf" )

#Create first PCA
x <- pcadapt(input = filename, K = 20)

#Plot the likelihoods
plot(x, option = "screeplot")
#Plot Plot the likelihoods for only first 10 K
plot(x, option = "screeplot", K = 10)

#Create population designations
poplist.names <- c(rep("POPA", 20),rep("POPB", 20),rep("POPC", 20), rep("POPD",20))

#Plot the actual PCA (first two PCAs)
plot(x, option = "scores", pop = poplist.names)
#Plot PCA with PCA 2 and PCA 3
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
#Plot PCA with PCA 3 and PCA 4
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

#Redo PCA with only 3 K
x <- pcadapt(filename, K = 3)

summary(x)

#Start looking for outliers
#Make Manhattan Plot
plot(x , option = "manhattan")
#Make qqplot
plot(x, option = "qqplot", threshold = 0.1)
# Look at P-value distribution
plot(x, option = "stat.distribution")

# Set FDR
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("qvalue")
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1

# Save outliers
outliers <- which(qval < alpha)


# Testing for library effects

poplist.names <- c(rep("LIB1", 40),rep("LIB2", 40))
x <- pcadapt(input = filename, K = 20)

plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)



x <- pcadapt(filename, K = 2)

summary(x)

plot(x , option = "manhattan")
plot(x, option = "qqplot", threshold = 0.1)

plot(x, option = "stat.distribution")

library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)

### Install required package devtools if not installed
if (!("devtools" %in% installed.packages())){install.packages("devtools")}
library(devtools)

### Install required package qvalue if not installed
if (!("qvalue" %in% installed.packages())){
  source("http://bioconductor.org/biocLite.R")
  biocLite("qvalue")
}

### Install OutFLANK if not installed
if (!("OutFLANK" %in% installed.packages())){install_github("whitlock/OutFLANK")}

if (!("vcfR" %in% installed.packages())){install.packages("vcfR")}
if (!("bigsnpr" %in% installed.packages())){install.packages("bigsnpr")}
library(OutFLANK)
library(vcfR)
library(bigsnpr)   # package for LD pruning

my_vcf <- read.vcfR("SNP.TRSdp5p05FHWE2A.recode.vcf")

geno <- extract.gt(my_vcf) # Character matrix containing the genotypes
position <- getPOS(my_vcf) # Positions in bp
chromosome <- getCHROM(my_vcf) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

G[is.na(G)] <- 9

head(G[,1:10])

pop <- read.table("popmap", header=FALSE)
pop <- pop$V2


my_fst <- MakeDiploidFSTMat(t(G), locusNames = paste0(chromosome,"_", position), popNames = pop)

my_dist <- OutFLANK(my_fst, NumberOfSamples = 4, qthreshold=0.1, RightTrimFraction=0.1, LeftTrimFraction=0.1)


OutFLANKResultsPlotter(my_dist)

plot(my_dist$results$FST, col=as.numeric(as.factor(chromosome)))


my_dist$results[which(my_dist$results$OutlierFlag == TRUE),]   

table_bay <- read.table("bayenv.final",header=TRUE)
plot(table_bay$BF1)

table_bay[which(table_bay$BF1 > 100),]
