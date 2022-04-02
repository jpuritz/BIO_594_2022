## BES594 assignment 7

This exercise is looking at calling and filter SNPs

#### Start by creating a conda environment

```
conda create -n week7 ddocent 
conda activate week7
```

Create working directory and link the data to it

```
mkdir week7
cd week7
ln -s /home/BIO594/Exercises/Week07_and_Week_08/simulated/*.fq.gz .
```

#### run [dDocent](https://www.ddocent.com)

```
dDocent 
```

Set parameters

```
80 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 80 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
yes
What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap.
Then press [ENTER]
PE
Reads will be assembled with Rainbow
CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa.
Would you like to enter a new c parameter now? Type yes or no and press [ENTER]
yes
Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)
0.9
Do you want to map reads?  Type yes or no and press [ENTER]
yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
4
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
6
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
jillashey@uri.edu


dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background.

Trimming reads and simultaneously assembling reference sequences

                                                                                                                        
                                                                                                                        
                       Number of Unique Sequences with More than X Coverage (Counted within individuals)                
   105000 +---------------------------------------------------------------------------------------------------------+   
          |           +           +          +           +           +           +          +           +           |   
          |*                                                                                                        |   
   100000 |-***                                                                                                   +-|   
          |    ****                                                                                                 |   
          |        ******                                                                                           |   
          |              ******                                                                                     |   
    95000 |-+                  ***********                                                                        +-|   
          |                               ******                                                                    |   
          |                                     ******                                                              |   
    90000 |-+                                         ******                                                      +-|   
          |                                                 ******                                                  |   
          |                                                       ******                                            |   
    85000 |-+                                                           ******                                    +-|   
          |                                                                   ******                                |   
          |                                                                         ******                          |   
    80000 |-+                                                                             *****                   +-|   
          |                                                                                    ******               |   
          |                                                                                          ******         |   
          |                                                                                                ******   |   
    75000 |-+                                                                                                    ***|   
          |                                                                                                         |   
          |           +           +          +           +           +           +          +           +           |   
    70000 +---------------------------------------------------------------------------------------------------------+   
          2           4           6          8           10          12          14         16          18          20  
                                                           Coverage                                                     
                                                                                                                        
Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly
5

                                                                                                                        
                                                                                                                        
                                 Number of Unique Sequences present in more than X Individuals                          
     3500 +---------------------------------------------------------------------------------------------------------+   
          |            +             +            +            +            +             +            +            |   
          |                                                                                                         |   
          |                                                                                                         |   
          |    *                                                                                                    |   
     3000 |-+   *                                                                                                 +-|   
          |      *                                                                                                  |   
          |       *                                                                                                 |   
          |        *                                                                                                |   
     2500 |-+       *                                                                                             +-|   
          |          *                                                                                              |   
          |           **                                                                                            |   
          |             **                                                                                          |   
          |               **                                                                                        |   
     2000 |-+               *****                                                                                 +-|   
          |                      ***                                                                                |   
          |                         *****                                                                           |   
          |                              ********                                                                   |   
     1500 |-+                                    ********                                                         +-|   
          |                                              *************                                              |   
          |                                                           *************                                 |   
          |                                                                        *******************              |   
          |            +             +            +            +            +             +           **************|   
     1000 +---------------------------------------------------------------------------------------------------------+   
          0            5             10           15           20           25            30           35           40  
                                                     Number of Individuals                                              
                                                                                                                        
Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals
3
```

Now dDocent will give this message:

```
At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

# As its running: 
dDocent assembled 2691 sequences (after cutoffs) into 1000 contigs

Using BWA to map reads

Creating alignment intervals
Genomic interval creation completed successfully.

Using FreeBayes to call SNPs
popmap file not found, building a new one from the sample names
100% 41:0=0s 14                                                                                                                                                                                     

Using VCFtools to parse TotalRawSNPS.vcf for SNPs that are called in at least 90% of individuals
```

okay well my analysis took 4 mins to run so I didn't get to sit back and relax. Here's the output

```
dDocent has finished with an analysis in /home/jashey/week7 

dDocent started Tue Mar 22 21:15:28 EDT 2022 

dDocent finished Tue Mar 22 21:19:33 EDT 2022 

After filtering, kept 3113 out of a possible 3268 Sites 
```

So we got 3113 SNPs to start. 

#### Filter some SNPs

Make a directory for the filtered data 

```
mkdir Filter
cd Filter/
ln -s ../TotalRawSNPs.vcf .
``` 

Use [vcftools](https://vcftools.github.io/index.html) to begin filtering. 

```
vcftools --vcf TotalRawSNPs.vcf --max-missing 0.5 --maf 0.001 --minQ 20 --recode --recode-INFO-all --out TRS
```

The parameters mean: 

Here's the output 

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TotalRawSNPs.vcf
	--recode-INFO-all
	--maf 0.001
	--minQ 20
	--max-missing 0.5
	--out TRS
	--recode
	
## gave some warnings to me, but it seems all good 

After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 3140 out of a possible 3268 Sites
Run Time = 1.00 seconds
```

So now we have 3140 SNPs. Now we are going to do more stringent filtering by lowering the end of the depth threshold. 

```
vcftools --vcf TRS.recode.vcf --minDP 5 --recode --recode-INFO-all --out TRSdp5
```

Here's the output. 

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRS.recode.vcf
	--recode-INFO-all
	--minDP 5
	--out TRSdp5
	--recode

After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 3140 out of a possible 3140 Sites
Run Time = 0.00 seconds
```

Still have 3140 SNPs. This makes sense because its simulated data! 

Okay now let's get rid of individuals that didn't sequence well. 

```
pop_missing_filter.sh TRSdp5.recode.vcf ../popmap 0.05 1 TRSdp5p05
```

Here's the output once again: 

```
After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 2946 out of a possible 3140 Sites
Run Time = 0.00 seconds
```

Now we are down to 2946 SNPs!

Let's use some dDocent filters now (as described in the dDocent filtering [tutorial](https://www.ddocent.com/filtering/)

```
dDocent_filters TRSdp5p05.recode.vcf TRSdp5p05
```

Here is the output

```
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 474 of 2946 

Are reads expected to overlap?  In other words, is fragment size less than 2X the read length?  Enter yes or no.
no
Number of additional sites filtered based on overlapping forward and reverse reads
 0 of 2472 

Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of additional sites filtered based on properly paired status
 0 of 2472 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 292 of 2472 


                                                                                                                        
                                                                                                                        
                                               Histogram of mean depth per site                                         
      300 +---------------------------------------------------------------------------------------------------------+   
          |   +    +    +     +    +     +    +    +     +    +    +     +    +     +    +    +     +    +     +    |   
          |                                               'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |   
          |                                                                                     *****               |   
      250 |-+                                                                                   * * *             +-|   
          |                                                                                     * * ***             |   
          |                                                                                     * * * *             |   
          |                                                                                  **** * * *             |   
      200 |-+                                                                                ** * * * *           +-|   
          |                                                                                  ** * * * **            |   
          |                                                                                **** * * * **            |   
      150 |-+                                                                            *** ** * * * **          +-|   
          |                                                                              * * ** * * * ****          |   
          |                                                                              * * ** * * * ** *          |   
          |                                                                              * * ** * * * ** *          |   
      100 |-+                                                                            * * ** * * * ** *        +-|   
          |                                                                              * * ** * * * ** *          |   
          |                                                                          *** * * ** * * * ** *          |   
          |                                                                          * *** * ** * * * ** *          |   
       50 |-+                                                                       ** * * * ** * * * ** *        +-|   
          |                                                                       **** * * * ** * * * ** *****      |   
          |                                                                       * ** * * * ** * * * ** * * ****   |   
          |   +    +    +     +    +     +    +    +     +    +    +     +    + *** ** * * * ** * * * ** * * * **** |   
        0 +---------------------------------------------------------------------------------------------------------+   
              12   15   18    21   24    27   30   33    36   39   42    45   48    51   54   57    60   63    66   69  
                                                          Mean Depth                                                    
                                                                                                                        
If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 5156.132535
The 95% cutoff would be 63
Would you like to use a different maximum mean depth cutoff than 63, yes or no
yes
Please enter new cutoff
75
Number of sites filtered based on maximum mean depth
 0 of 2472 

Number of sites filtered based on within locus depth mismatch
 0 of 2250 

Total number of sites filtered
 696 of 2946 

Remaining sites
 2250 

Filtered VCF file is called TRSdp5p05.FIL.recode.vcf

Filter stats stored in TRSdp5p05.filterstats
```

So it looks like we have 2250 remaining SNPs??

We now change the formatting of vcf file and save as a prim file. Then we can run the prim file through the vcftools software to remove the insertions/deletions

```
vcfallelicprimitives -k -g TRSdp5p05.FIL.recode.vcf |sed 's:\.|\.:\.\/\.:g' > TRSdp5p05F.prim
vcftools --vcf TRSdp5p05F.prim --recode --recode-INFO-all --remove-indels --out SNP.TRSdp5p05F
```

Output is here: 

```
After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 1888 out of a possible 2431 Sites
Run Time = 0.00 seconds
```

Now we got 1888 SNPs.

Next we filter for the H-W Equilibrium. Lots of filtering 

```
filter_hwe_by_pop.pl -v SNP.TRSdp5p05F.recode.vcf -p ../popmap -c 0.5 -out SNP.TRSdp5p05FHWE
```

Here's the output: 

```
Processing population: PopA (20 inds)
Processing population: PopB (20 inds)
Processing population: PopC (20 inds)
Processing population: PopD (20 inds)
Outputting results of HWE test for filtered loci to 'filtered.hwe'
Kept 1888 of a possible 1888 loci (filtered 0 loci)
```

Kept all SNPs, makes sense because this is simulared data 

Now filter to keep SNPs with minor allele frequencies of 0.05 and higher

```
After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 925 out of a possible 1888 Sites
Run Time = 0.00 seconds
```

Output for that: 

```
After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 925 out of a possible 1888 Sites
Run Time = 0.00 seconds
```

We are now at 925 SNPs!

Let's convert files from VCF to other formats using [PGDspider](http://www.cmpg.unibe.ch/software/PGDSpider/) b/c other analyses that will be done later need certain file formats. 

Copy the PGDspider configuration file and file to map individuals to population

```
cp /home/BIO594/DATA/Week7/example/BSsnp.spid .
ln -s ../popmap .
```

Run PGDspider (which uses java)

```
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRSdp5p05FHWEmaf05.recode.vcf -outputfile SNP.TRSdp5p05FHWEBS -spid BSsnp.spid
```

Here's the output:

```
initialize convert process...
read input file...
read input file done.
write output file...
write output file done.
```

#### Outlier detection 

[BayeScan](http://cmpg.unibe.ch/software/BayeScan/index.html) attempts to find loci under selection using differences in allele frequencies between populations. Putative outliers are measured with Fst coefficients.

In this code, nbp is the "number of pilot runs" (default 20) and thin is the "thinning interval size" or the number of iterations between two samples (default 10). These values set the parameters for the Markov "chain". Because this code is an interative process, it takes a while to run.

Run BayeScan 

```
BayeScan2.1_linux64bits SNP.TRSdp5p05FHWEBS -nbp 30 -thin 20
```

Hmm it says its using 80 threads. Hopefully not using all the KITT RAM rn.

Output from BayeScan run

```
Using 80 threads (80 cpu detected on this machine)
Pilot runs...
1% 2% 3% 4% 5% 6% 7% 8% 9% 10% 11% 12% 13% 14% 15% 16% 17% 18% 19% 20% 21% 22% 23% 24% 25% 26% 27% 28% 29% 30% 31% 32% 33% 34% 35% 36% 37% 38% 39% 40% 41% 42% 43% 44% 45% 46% 47% 48% 49% 50% 51% 52% 53% 54% 55% 56% 57% 58% 59% 60% 61% 62% 63% 64% 65% 66% 67% 68% 69% 70% 71% 72% 73% 74% 75% 76% 77% 78% 79% 80% 81% 82% 83% 84% 85% 86% 87% 88% 89% 90% 91% 92% 93% 94% 95% 96% 97% 98% 99% 
Calculation...
1% 2% 3% 4% 5% 6% 7% 8% 9% 10% 11% 12% 13% 14% 15% 16% 17% 18% 19% 20% 21% 22% 23% 24% 25% 26% 27% 28% 29% 30% 31% 32% 33% 34% 35% 36% 37% 38% 39% 40% 41% 42% 43% 44% 45% 46% 47% 48% 49% 50% 51% 52% 53% 54% 55% 56% 57% 58% 59% 60% 61% 62% 63% 64% 65% 66% 67% 68% 69% 70% 71% 72% 73% 74% 75% 76% 77% 78% 79% 80% 81% 82% 83% 84% 85% 86% 87% 88% 89% 90% 91% 92% 93% 94% 95% 96% 97% 98% 99% 
```

To look at the plot, we need to copy the R source code.

```
cp /home/BIO594/DATA/Week7/example/plot_R.r .
```

Okay now we use R. I totally forget how to use R in terminal/ on KITT. also don't remember VSC method. So maybe I'll just copy the files to my computer and do it in my own RStudio???

My brain is broken, ask Megan tomorrow. 

just copied it onto my own laptop for R. Ran this code in R 

```
source("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plot_R.r")
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plot/BS_plot.png") 
plot_bayescan("~/Desktop/URI/Spring2022/BIO594/exercises/week7/SNP.TRSdp5p05FH_fst.txt")

$outliers
[1] 799 800 874 875 895

$nb_outliers
[1] 5
```

Which generated this plot: 

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/BS_plot.png)

Now we need to limit SNPs to only those with two alleles with this code: 

```
vcftools --vcf SNP.TRSdp5p05FHWEmaf05.recode.vcf --max-alleles 2 --recode --recode-INFO-all --out SNP.TRSdp5p05FHWE2A
```

Output here: 

```
After filtering, kept 80 out of 80 Individuals
Outputting VCF file...
After filtering, kept 923 out of a possible 925 Sites
Run Time = 0.00 seconds
```

Now we work with [PCAdapt](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html) in R. Again just copied files from KITT server to my computer 

```{r}
#Load pcadapt library
library(pcadapt)
library(vcfR)

#load our VCF file into R
filename <- read.pcadapt("~/Desktop/URI/Spring2022/BIO594/exercises/week7/SNP.TRSdp5p05FHWE2A.recode.vcf", type = "vcf" )

#Create first PCA
x <- pcadapt(input = filename, K = 20)

#Plot the likelihoods
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/likelihoods.png")
plot(x, option = "screeplot")
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/likelihoods.png)

```
#Plot the likelihoods for only first 10 K
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/likelihoods_10K.png")
plot(x, option = "screeplot", K = 10)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/likelihoods_10K.png)

```
#Create population designations
poplist.names <- c(rep("POPA", 20),rep("POPB", 20),rep("POPC", 20), rep("POPD",20))

#Plot the actual PCA (first two PCAs)
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/PC1_PC2.png")
plot(x, option = "scores", pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/PC1_PC2.png)

```
#Plot PCA with PCA 2 and PCA 3
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/PC2_PC3.png")
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/PC2_PC3.png)

```
#Plot PCA with PCA 3 and PCA 4
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/PC3_PC4.png")
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/PC3_PC4.png)

```
#Redo PCA with only 3 K
x <- pcadapt(filename, K = 3)
summary(x)

#Plot the actual PCA (first two PCAs)
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/3K_PC1_PC2.png")
plot(x, option = "scores", pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/3K_PC1_PC2.png)

```
#Plot PCA with PCA 2 and PCA 3
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/3K_PC2_PC3.png")
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/3K_PC2_PC3.png)

```
#Start looking for outliers
#Make Manhattan Plot
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/manhattan.png")
plot(x , option = "manhattan")
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/manhattan.png)

```
#Make qqplot
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/qqplot.png")
plot(x, option = "qqplot", threshold = 0.1)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/qqplot.png)

```
# Look at P-value distribution
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/pval.png")
plot(x, option = "stat.distribution")
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/pval.png)

```
# Set FDR
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1

# Save outliers
outliers <- which(qval < alpha)

# Testing for library effects

poplist.names <- c(rep("LIB1", 40),rep("LIB2", 40))
x <- pcadapt(input = filename, K = 20)
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/library-PC1_PC2.png")
plot(x, option = "scores", pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/library-PC1_PC2.png)

```
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/library-PC2_PC3.png")
plot(x, option = "scores", i = 2, j = 3, pop = poplist.names)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/library-PC2_PC3.png)

```
x <- pcadapt(filename, K = 2)

summary(x)
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/2K-manhattan.png")
plot(x , option = "manhattan")
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/2K-manhattan.png)

```
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/2K-qqplot.png")
plot(x, option = "qqplot", threshold = 0.1)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/2K-qqplot.png)

```
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/2K-pval.png")
plot(x, option = "stat.distribution")
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/2K-pval.png)

```
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
```

Run [outflank](https://github.com/whitlock/OutFLANK) in R

```
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)   # package for LD pruning

my_vcf <- read.vcfR("~/Desktop/URI/Spring2022/BIO594/exercises/week7/SNP.TRSdp5p05FHWE2A.recode.vcf")

# Scanning file to determine attributes.
# File attributes:
#   meta lines: 65
#   header_line: 66
#   variant count: 923
#   column count: 89
# Meta line 65 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
#   Character matrix gt rows: 923
#   Character matrix gt cols: 89
#   skip: 0
#   nrows: 923
#   row_num: 0
# Processed variant: 923
# All variants processed

geno <- extract.gt(my_vcf) # Character matrix containing the genotypes
position <- getPOS(my_vcf) # Positions in bp
chromosome <- getCHROM(my_vcf) # Chromosome information

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

G[is.na(G)] <- 9

head(G[,1:10])

#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    2    2    1    2    2    2    2    2    2     2
# [2,]    2    1    2    1    1    2    2    1    2     2
# [3,]    2    2    1    1    1    2    2    2    2     2
# [4,]    0    0    0    0    0    1    1    0    0     0
# [5,]    0    0    0    0    0    0    0    0    0     0
# [6,]    0    0    0    0    0    0    0    0    0     0

pop <- read.table("~/Desktop/URI/Spring2022/BIO594/exercises/week7/popmap", header=FALSE)
pop <- pop$V2

# calculates Fst coefficients
my_fst <- MakeDiploidFSTMat(t(G), locusNames = paste0(chromosome,"_", position), popNames = pop)

my_dist <- OutFLANK(my_fst, NumberOfSamples = 4, qthreshold=0.1, RightTrimFraction=0.1, LeftTrimFraction=0.1)
```

```
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/outflank-fst_without_sample_size_correction.png")
OutFLANKResultsPlotter(my_dist)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/outflank-fst_without_sample_size_correction.png)

```
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/csome_plot.png")
plot(my_dist$results$FST, col=as.numeric(as.factor(chromosome)))
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/csome_plot.png)

```
my_dist$results[which(my_dist$results$OutlierFlag == TRUE),]
```

Now go back to terminal to run some [BayEnv2](https://github.com/lotteanna/Bioinformatics/blob/master/GBS_BayEnv.md). convert vcf file into BayEnv output 

```
cp /home/BIO594/DATA/Week7/example/SNPBayEnv.spid .
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.TRSdp5p05FHWE2A.recode.vcf -outputfile SNP.TRSdp5p05FHWEBayEnv.txt -spid SNPBayEnv.spid
```

Output: 

```
WARN  14:53:26 - PGDSpider configuration file not found! Loading default configuration.
initialize convert process...
read input file...
read input file done.
write output file...
write output file done.
```

now run BayEnv to generate covariance matrix 

```
bayenv2 -i SNP.TRSdp5p05FHWEBayEnv.txt -p 4 -k 100000 -r 63479 > matrix.out
```

pull out only the last iteration

```
tail -5 matrix.out | head -4 > matrix
```

With the matrix we will use our environmental factor file:

```
cat environ

-0.888330138	-0.565300997	0.080757285	1.37287385
-0.565300997	-0.484543712	-0.565300997	-0.403786427
```

You have to make this file! This is your environmental data split up by population. For this example, I used the following data. Be sure that your values are separated by tabs!

The environmental file are standardized environmental data with each line representing an environemtal factor with the value for each population tab delimited. The dummy file has 2 variables for 4 populations

Now calculate the Bayes Factor for each SNP for each environmental variable:

```
ln -s /usr/local/bin/bayenv2 .
calc_bf.sh SNP.TRSdp5p05FHWEBayEnv.txt environ matrix 4 10000 2
```

took a few mins. Now convert bayenv to format that's good for R

```
paste <(seq 1 923) <(cut -f2,3 bf_environ.environ ) > bayenv.out
cat <(echo -e "Locus\tBF1\tBF2") bayenv.out > bayenv.final
```

Go back to R and plot bayenv 

```
table_bay <- read.table("~/Desktop/URI/Spring2022/BIO594/exercises/week7/bayenv.final",header=TRUE)
png("~/Desktop/URI/Spring2022/BIO594/exercises/week7/plots/bayenv2_plot.png")
plot(table_bay$BF1)
dev.off()
```

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/Week07/Jill/plots/bayenv2_plot.png)

```
table_bay[which(table_bay$BF1 > 100),]
```

Output here: 
```
    Locus    BF1     BF2
873   873 201.66 0.66737
```
