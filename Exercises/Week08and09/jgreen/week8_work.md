# Week 8 and 9 problems sets

## Assignment breakdown

## Workflow

### Part 1

* Import data and scripts into a working directory
* Run BayeScan on each dataset
* Plot Bayescan in Rstudio
* Refilter based on distribution
* Intial plotting of SNPs in Rstudio using PCadapt
* Filter outliers using qvalue package
* Replot SNPs
* Apply OutFLANK to the dataset and plot Fst
* Run BayEnv2 on final filter set
* Conduct a PCA analysis on outlier free dataset
* Conduct a DAPC on an outlier dataset
* Perform two analysis from the Silliman et al 2019

### Part 2

* Repeat the Sillman Analysis
* Repeat the Zyck Analysis


## Part 1 Starting information

In the directory, `/home/BIO594/Exercises/Week07_and_Week_08`, you will find a few useful files and two directories, `realdata` and `simulated`

I would like you to complete the following and push it to your own folder in this repository directory.

* For the realdata
  * Run BayeScan
  * Run BayEnv
  * Run at least one PCA and one DAPC using outlier free data set
  * Perform at least two analyses from Silliman et al	
    
Each directory has a `popmap` `environ` and `LibraryInfo` files

`popmap` is a mapping of individuals to populations

`environ` is an environmental factor file for BayEnv

`LibraryInfo` maps samples to sequencing library

There are also PGDspider configuration files and an R file for BayeScan output plotting.

## Part 2

* Repeat the redundancy analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/RedundancyAnalysis/RDA_Outlier_Hap.Rmd)
* Repeat the EEMS analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/EEMS/NB_EEMS_OutlierHap.md)

## Workflow

### Part 1.1

### Realdata set

* Run BayeScan on the realdataset
* Plot Bayescan in Rstudio
* Refilter based on distribution
* Intial plotting of SNPs in Rstudio using PCadapt
* Filter outliers using qvalue package
* Replot SNPs
* Apply OutFLANK to the dataset and plot Fst
* Run BayEnv2 on final filter set
* Conduct a PCA analysis on outlier free dataset
* Conduct a DAPC on an outlier dataset

### Part 1.2

* Repeat the Sillman Analysis
* Repeat the Zyck Analysis

### Part 2

* Repeat the redundancy analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/RedundancyAnalysis/RDA_Outlier_Hap.Rmd)
* Repeat the EEMS analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/EEMS/NB_EEMS_OutlierHap.md)

### Setting up the working environment

Make a directory for week 8 work

```
cd ~/repos/BIO594_work/
mkdir Week8
```

Create a conda env for Week 8 work uploading ddocent

```
conda create -n week8 ddocent
conda activate week8
```

Link the directories for each `realdata` and `simulated` directories

```
cd Week8
cp -r /home/BIO594/Exercises/Week07_and_Week_08/realdata/ .
```
### Part 1 Realdataset

Convert the recoded vcf file using pDGspider

```
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf -outputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL -spid BSsnp.spid
```

Use Bayescan to detect Fst outliers within the dataset

```
BayeScan2.1_linux64bits SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL -threads 20 -nbp 30 -thin 20
```

```R
source("plot_R.r")
plot_bayescan("SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE_fst.txt")
```

![](realdata/figures/initial_bayescan.png)

```R
#Load pcadapt library
install.packages("pcadapt", dependencies = TRUE)
library(pcadapt)

#load our VCF file into R
filename <- read.pcadapt("SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf", type = "vcf" )
```

1816 variant(s) have been discarded as they are not SNPs.
Summary:

	- input file:				SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf
	- output file:				/tmp/RtmpVkjmxP/file8b4833b5a309.pcadapt

	- number of individuals detected:	384
	- number of loci detected:		7303

5487 lines detected.
384 columns detected.

```R
#Create first PCA
x <- pcadapt(input = filename, K = 20)

#Plot the likelihoods
plot(x, option = "screeplot")
```

![](realdata/figures/k20_screeplot.png)

```R
#Plot Plot the likelihoods for only first 10 K
plot(x, option = "screeplot", K = 10)
```

![](realdata/figures/k10_screeplot.png)

```R
#Create population designations
poplist.names <- c(rep("ACM", 19),rep("BHS", 19),rep("BRS", 26), rep("CCM",70), rep("DMS",16), rep("DRM",11), rep("FMS",19), rep("IPM",12), rep("LHM",62), rep("MCM",11), rep("PGS",20), rep("PMM",20), rep("SPS",20), rep("SSM",20), rep("TCS",20), rep("WPS",19))

#Plot the actual PCA (first two PCAs)
plot(x, option = "scores", pop = poplist.names)
```

![](realdata/figures/pcc1_pc2.png)

```R
#Plot PCA with PCA 2 and PCA 3
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
```

![](realdata/figures/pc2_pc3.png)

```R
#Plot PCA with PCA 3 and PCA 4
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
```

![](realdata/figures/pc3_pc4.png)

```R
#Redo PCA with only 3 K
x <- pcadapt(filename, K = 3)

summary(x)

                Length Class  Mode   
scores           1152  -none- numeric
singular.values     3  -none- numeric
loadings        16461  -none- numeric
zscores         16461  -none- numeric
af               5487  -none- numeric
maf              5487  -none- numeric
chi2.stat        5487  -none- numeric
stat             5487  -none- numeric
gif                 1  -none- numeric
pvalues          5487  -none- numeric
pass             3092  -none- numeric
```

```R
#Start looking for outliers
#Make Manhattan Plot
plot(x , option = "manhattan")
```

![](realdata/figures/manhattan_k20.png)

```R
#Make qqplot
plot(x, option = "qqplot", threshold = 0.1)
```

![](realdata/figures/qqplot_k20.png)

```R
# Look at P-value distribution
plot(x, option = "stat.distribution")
```

![](realdata/figures/stat-dist_k20.png)


```R
# Set FDR
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1

# Save outliers
outliers <- which(qval < alpha)
invisible(lapply(outliers, write, "outliers.txt", append=TRUE))
```

Count how many outliers there are in the .txt file

```text
cat outliers.txt | wc - w
99
```

Save outlier to new file

```
mawk '!/#/' SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf | cut -f1,2 > totalloci

NUM=(`cat totalloci | wc -l`)

paste <(seq 1 $NUM) totalloci > loci.plus.index

cat outliers.txt | parallel "grep -w ^{} loci.plus.index" | cut -f2,3> outlier.loci.txt

head outliter.loci.txt
```
# Testing for library effects
poplist.names <- c(rep("LIB1", 40),rep("LIB2", 40))
x <- pcadapt(input = filename, K = 20)

plot(x, option = "scores", pop = poplist.names)
```

![](figure/outlier_pc1_pc2.png)


```R
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
```

![](figure/outliers_pc2_pc3.png)

```R
x <- pcadapt(filename, K = 2)

summary(x)

                Length Class  Mode   
scores           160   -none- numeric
singular.values    2   -none- numeric
loadings        1530   -none- numeric
zscores         1530   -none- numeric
af               765   -none- numeric
maf              765   -none- numeric
chi2.stat        765   -none- numeric
stat             765   -none- numeric
gif                1   -none- numeric
pvalues          765   -none- numeric
pass             765   -none- numeric
```

```R
# Set FDR
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("qvalue")
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)

### Install required package devtools if not installed
#if (!("devtools" %in% installed.packages())){install.packages("devtools")}
library(devtools)

### Install required package qvalue if not installed
#if (!("qvalue" %in% installed.packages())){ source("http://bioconductor.org/biocLite.R") biocLite("qvalue") }

### Install OutFLANK if not installed
#if (!("OutFLANK" %in% installed.packages())){install_github("whitlock/OutFLANK")}

#if (!("vcfR" %in% installed.packages())){install.packages("vcfR")}
#if (!("bigsnpr" %in% installed.packages())){install.packages("bigsnpr")}
library(OutFLANK)
library(vcfR)
library(bigsnpr)   # package for LD pruning

my_vcf <- read.vcfR("SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf")

Scanning file to determine attributes.
File attributes:
  meta lines: 64
  header_line: 65
  variant count: 7303
  column count: 393
Meta line 64 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 7303
  Character matrix gt cols: 393
  skip: 0
  nrows: 7303
  row_num: 0
Processed variant: 7303
All variants processed

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
```

![](realdata/figures/fst_wo_sample_size_correctio.png)

```R
plot(my_dist$results$FST, col=as.numeric(as.factor(chromosome)))
```

![](realdata/figures/dist_fst.png)

```R
my_dist$results[which(my_dist$results$OutlierFlag == TRUE),] 

LocusName        He          FST           T1         T2  FSTNoCorr    T1NoCorr   T2NoCorr meanAlleleFreq
31         Sc28pcJ_63_HRSCAF_70_732289 0.1152662  0.130247551  0.007666183 0.05885856 0.22556411 0.013339852 0.05913996      0.9385965
1731   Sc28pcJ_520_HRSCAF_571_10014978 0.3705361  0.042170872  0.007861199 0.18641302 0.08118297 0.015150515 0.18662184      0.7544248
2346    Sc28pcJ_649_HRSCAF_705_5313506 0.1788224  0.043527322  0.003921960 0.09010341 0.10264889 0.009267559 0.09028407      0.9007353
2515   Sc28pcJ_649_HRSCAF_705_20192478 0.4695168  0.090081335  0.021511949 0.23880584 0.17948779 0.043106804 0.24016566      0.6234568
2678   Sc28pcJ_649_HRSCAF_705_29016814 0.4938272 -0.206557377 -0.054687500 0.26475694 0.47693119 0.130859375 0.27437789      0.4444444
4362 Sc28pcJ_1199_HRSCAF_1291_10833636 0.1723356  0.070861651  0.006327986 0.08930057 0.28203396 0.025980392 0.09211796      0.0952381
5665 Sc28pcJ_1650_HRSCAF_1766_13774262 0.4758666 -0.007028933 -0.001678514 0.23880067 0.06474452 0.015503736 0.23946019      0.6098485
6176  Sc28pcJ_1840_HRSCAF_2023_1808587 0.4534616 -0.040147671 -0.009152118 0.22796138 0.11016296 0.025235293 0.22907239      0.6525424
6315  Sc28pcJ_1840_HRSCAF_2023_9754856 0.4591837 -0.138613861 -0.033333333 0.24047619 0.18006431 0.044444444 0.24682540      0.3571429
6889  Sc28pcJ_1843_HRSCAF_2042_6801215 0.2592670  0.019225304  0.002522507 0.13120763 0.17194876 0.022758078 0.13235384      0.8469388
     indexOrder GoodH      qvalues      pvalues pvaluesRightTail OutlierFlag
31           31 goodH 0.000000e+00 0.000000e+00     0.000000e+00        TRUE
1731       1731 goodH 2.356936e-03 1.044947e-05     5.224736e-06        TRUE
2346       2346 goodH 1.353295e-05 5.333181e-08     2.666590e-08        TRUE
2515       2515 goodH 0.000000e+00 0.000000e+00     0.000000e+00        TRUE
2678       2678 goodH 0.000000e+00 0.000000e+00     0.000000e+00        TRUE
4362       4362 goodH 0.000000e+00 0.000000e+00     0.000000e+00        TRUE
5665       5665 goodH 9.444949e-02 4.652684e-04     2.326342e-04        TRUE
6176       6176 goodH 2.291889e-06 7.903067e-09     3.951534e-09        TRUE
6315       6315 goodH 0.000000e+00 0.000000e+00     0.000000e+00        TRUE
6889       6889 goodH 2.253753e-13 6.661338e-16     3.330669e-16        TRUE
```

### Run BayEnv

Create a file that is converted vcf to .txt file using PGDspider 

```
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf -outputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.BayEnv.txt -spid SNPBayEnv.spid
```

Run BayEnv to generate the covariance matrix

```
bayenv2 -i SNP.TRSdp5p05FHWEBayEnv.txt -p 16 -k 100000 -r 63479 > matrix.out
```

This code generates 100,000 iterations.  We only need the last one.

```
tail -5 matrix.out | head -4 > matrix
```

With the matrix we will use our environmental factor file:

```text
cat environ

-6.23396E-06    -0.001053539    -0.001053539    -0.000903924    0.001090943     -0.000704437    -0.001452512    -0.000654566    -0.001651999    -0.001103411    -6.23396E-06    -0.001053539    -0.000155849    0.0088210510.001090943     -0.001203154
0.030764586     0.017797952     0.022785119     0.044728653     0.001340301     0.013808218     0.024281269     0.012810785     0.030265869     0.007324901     -0.000654566    0.019792818     0.020291535     0.0073249010.003335168     0.021787685
```

The environmental file are standardized environmental data with each line representing an environemtal factor with the value for each population tab delimited.  This dummy file has 2 variables for 16 populations

Next, we calculate the Bayes Factor for each SNP for each environmental variable:

```
ln -s /usr/local/bin/bayenv2 .
calc_bf.sh SNP.TRSdp5p05FHWEBayEnv.txt environ matrix 4 10000 2
```

Next, we convert the output into something suitable to input into R

```text
paste <(seq 1 7303) <(cut -f2,3 bf_environ.environ ) > bayenv.out
cat <(echo -e "Locus\tBF1\tBF2") bayenv.out > bayenv.final
```

What is in this bayenv.final file?

```bash
cat bayenv.final 
```

Now, open R and make sure you transfer the bayenv.final to the proper working directory

```R
table_bay <- read.table("bayenv.final",header=TRUE)
plot(table_bay$BF1)

table_bay[which(table_bay$BF1 > 100),]
```

# Run PCA in R from Class

```R
library(adegenet)
library(vcfR)


my_vcf <- read.vcfR("SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf")
```

Scanning file to determine attributes.
File attributes:
  meta lines: 64
  header_line: 65
  variant count: 7303
  column count: 393
Meta line 64 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 7303
  Character matrix gt cols: 393
  skip: 0
  nrows: 7303
  row_num: 0
Processed variant: 7303
All variants processed

```R
my_genind <- vcfR2genind(my_vcf)

strata<- read.table("LibraryInfo", header=TRUE)

strata_df <- data.frame(strata)

strata(my_genind) <- strata_df

setPop(my_genind) <- ~Population
```

```R
#Test Population Structure
library("hierfstat")
## fstat(my_genind)
summary(my_genind)

names(my_genind)
```

 [1] "tab"       "loc.fac"   "loc.n.all" "all.names" "ploidy"    "type"     
 [7] "other"     "call"      "pop"       "strata"    "hierarchy"

 ```R
pop(my_genind)
```

Levels: ACM BHS BRS CCM DMS DRM FMS IPM LHM MCM PGS PMM SPS SSM TCS WPS

```R
basic.stats(my_genind)

$overall
    Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.2025 0.2042 0.2043 0.0001 0.2043 0.0001 0.0003 0.0003 0.0087 0.0001 
```

```R
X <- tab(my_genind, freq = TRUE, NA.method = "mean")

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
```

![](realdata/figures/pca_eigenvalues.png)


```R
s.class(pca1$li, pop(my_genind))

title("PCA of simulated dataset\naxes 1-2")

add.scatter.eig(pca1$eig[1:20], 3,1,2)
```

![](realdata/figures/pca_real_1_2.png)

```R
col <- funky(15)

s.class(pca1$li, pop(my_genind),xax=1,yax=2, col=col, axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
```

![](realdata/figures/pca_location.png)

DAPC

```R
grp <- find.clusters(my_genind)

table(pop(my_genind), grp$grp)

table.value(table(pop(my_genind), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("ori", 1:16))
```

![](realdata/figures/grp_table_16_10.png)

```R
dapc1 <- dapc(my_genind, grp$grp)
Choose the number PCs to retain (>=1): 
4
Choose the number discriminant functions to retain (>=1): 
4
```

![](realdata/figures/dis_analysis_eigen_1.png)

```R
scatter(dapc1,col=col,bg="white", solid=1)
```

![](realdata/figures/scatter_dapc.png)

```R
contrib <- loadingplot(dapc1$var.contr, axis=1, thres=.01, lab.jitter=1)
```

![](realdata/figures/loading_plot_contrib1.png)

```R
contrib

setPop(my_genind) <- ~Library

dapc1 <- dapc(my_genind, pop(my_genind))

contrib <- loadingplot(dapc1$var.contr, axis=1, thres=.05, lab.jitter=1)
```

![](realdata/figures/loading_plot_2.png)

Structure Like

```R
compoplot(dapc1, posi="bottomright",txt.leg=paste("Cluster", 1:2), lab="", ncol=1, xlab="individuals")
```

![](realdata/figures/cluster_graph.png)

```R
temp <- which(apply(dapc1$posterior,2, function(e) all(e<0.9)))

compoplot(dapc1, subset=temp, posi="bottomright", txt.leg=paste("Cluster", 1:2), ncol=2)
```

![](realdata/figures/cluster_subtemp.png)

## Part 2