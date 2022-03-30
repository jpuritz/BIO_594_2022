## Assignment description:

### Part 1
In the directory, /home/BIO594/Exercises/Week07_and_Week_08, you will find a few useful files and two directories, `realdata` and `simulated`. I would like you to complete the following and push it to your own folder in this repository directory.

For the `realdata`
- Run BayeScan
- Run BayEnv
- Run at least one PCA and one DAPC using outlier free data set
- Perform at least two analyses from Silliman et al
*Each directory has a `popmap` `environ` and `LibraryInfo` files*

`popmap` is a mapping of individuals to populations

`environ` is an environmental factor file for BayEnv

`LibraryInfo` maps samples to sequencing library

There are also PGDspider configuration files and an R file for BayeScan output plotting.

### Part 2
Repeat the redundancy analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/RedundancyAnalysis/RDA_Outlier_Hap.Rmd)
Repeat the EEMS analysis documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/EEMS/NB_EEMS_OutlierHap.md)

# Beginning analysis on `real data`
For this exercise we are starting from the filtered VCF file. 

Create a new working directory then create a link to the `realdata` we will use in this analysis.
```
mkdir week8_9
cd week8_9
cp /home/BIO594/Exercises/Week07_and_Week_08/realdata/* .
```

## Converting from VCF to other outputs
Before moving on to our analyses, we will have to first convert from VCF to other file formats using PGDspider.

To do so, you will need a PGDspider configuration file and file to map individuals to population. For this instance, both the `BSsnp.spid` and the `popmap` files where provided.

Now, run PGDspider. This program is written in java and converts the vcf into a PDG file format that can be run in BayeScan.
The input file is your VCF.
```
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf -outputfile SNP.TRSdp5p05FHWEBS -spid BSsnp.spid
```

## Bayescan
### Run BayeScan
In brief, [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) attempts to find loci under selection using differences in allele frequencies between populations. Putative outliers are measured with Fst coefficients.

In this code, `nbp` is the "number of pilot runs" (default 20) and `thin` is the "thinning interval size" or the number of iterations between two samples (default 10). These values set the parameters for the Markov "chain". Because this code is an interative process, it takes a *very long time* to run ()
```
BayeScan2.1_linux64bits SNP.TRSdp5p05FHWEBS -nbp 30 -thin 20
```

Because this takes so long, you can keep it running in the background using `Cmd+Z` then type `bg`, `ENTER`. 

To disown, type `disown -a`. 

To check the status of what is running in the background, you can use `top` to look at everything running on kitt. To check out what you're specifically running, you can use `jobs`.

After disowning, you will not see your job under `jobs` list but it will still be running on kitt when you check out `top`.

If you know you're going to be running a job for a while, you could also use `tmux`. 

Make sure you have the R source file `plot_R.r` in your current working directory. In this example, it was provided in the `realdata` folder. 

### Open Rstudio and plot BayeScan output. -- next step as of 3/30/22 @ 5:30pm
Run the following to plot the output of BayeScan.
```
source("plot_R.r")
png("BS_plot.png") 
plot_bayescan("SNP.TRSdp5p05FH_fst.txt")
dev.off()
```

Resulting plot
![](figures/BS_plot.png)


## BayEnv2
[Documentation](https://bitbucket.org/tguenther/bayenv2_public/src)

This analysis will start on the terminal (in the week8_9 directory) and will transition to R studio at the end.

First, you have to convert vcf to BayEnv input with the following code.
```
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile SNP.DP3g98maf01_85INDoutFIL.NO2a.HWE.FIL.recode.vcf -outputfile SNP.TRSdp5p05FHWEBayEnv.txt -spid SNPBayEnv.spid
```

Run BayEnv to generate the covariance matrix

`bayenv2 -i SNP.TRSdp5p05FHWEBayEnv.txt -p 4 -k 100000 -r 63479 > matrix.out`

This code generates 100,000 iterations.  We only need the last one.

`tail -5 matrix.out | head -4 > matrix`

With the matrix we will use our environmental factor file:

`cat environ` 

**You have to make this file! This is your environmental data split up by population. For this example, the file was provided in the `realdata` directory. Be sure that your values are separated by *tabs*!**

The environmental file are standardized environmental data with each line representing an environemtal factor with the value for each population tab delimited.  

Next, we calculate the Bayes Factor for each SNP for each environmental variable:

```
ln -s /usr/local/bin/bayenv2 .
calc_bf.sh SNP.TRSdp5p05FHWEBayEnv.txt environ matrix 4 10000 2
```

Next, we convert the output into something suitable to input into R
```
paste <(seq 1 923) <(cut -f2,3 bf_environ.environ ) > bayenv.out
cat <(echo -e "Locus\tBF1\tBF2") bayenv.out > bayenv.final
```

Now, you can go into R studio and run the following...

`R`

```
table_bay <- read.table("bayenv.final",header=TRUE)
png("bayenv2_plot.png")
plot(table_bay$BF1)
dev.off()
```
![](figures/bayenv2_plot.png)

```
table_bay[which(table_bay$BF1 > 100),]
```
Output:

```
    Locus    BF1     BF2
872   872 508.46 3.47660
873   873 603.51 0.51944
```