# Analyzing Pumpkin's DNA
## WGS Data Analysis Pipeline
![Pumpkin](https://github.com/jpuritz/BIO_594_2022/tree/main/Exercises/course_project/mgregoire/pumpkin_resize.jpg)
Pumpkin got his DNA sequenced through Basepaws! Basepaws delivered a breed profile, dental report, and genetic report of Pumpkin's SNPs. Below is a pipeline I have established for some analyses that can be performed on the resulting fastq and vcf files provided by Basepaws.

## Set up directories, softwares, and files 
Make sure that bioconda is is installed. I am using the Kitt server and bioconda has already been installed. Bioconda makes it easy to install many of the various programs that will be used in this pipeline (eg: fastqc).
Make a folder for the project and cd into it. 
- `mkdir FinalProject | cd FinalProject`

Since the files were supplied to me on a flashdrive from Basepaws, I sftp'd the files onto the kitt server with the following commands:
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd D: | lls`
- `put AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq.gz | put AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq.gz | put AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq.gz | put pumpkin_102.hard-filtered.gvcf | put pumpkin_102.hard-filtered.gvcf.gz.tbi`

I wish I knew what the information in these file names mean! Since I don't know, we'll go ahead and analyze all 3 files.

Unzip the files
- `gunzip AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq.gz`
- `gunzip AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq.gz`
- `gunzip AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq.gz`

## Run fastqc on the files to check for quality
Run fastqc and download the .html output files
- `fastqc *.fastq`
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Documents\Pumpkin`
-  `get AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS_fastqc.html | get AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS_fastqc.html | get AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS_fastqc.html`

The fastqc files from Pumpkin's DNA reveal the following:
- **Per base sequence quality**: The quality scores across all bases look good for all 3 fastq files! The reads start to drop off towards the tail ends (this is expected as the reads get longer), but the scores stay in the green zone the entire time hovering around 36-34 and only with error bars dropping to about 25 from 100-150bp. This means that the average of the reads had >99.9% accuracy.
- **Per sequence quality scores**: The quality score distribution across all sequences peak at 36 for all three files. This means that Q36 had more read numbers than other quality scores.
- **Per base sequence content**: This metric shows that the percentage of each of the four nucleotides (A, T, C, G) at each position across all reads showed normal random distributions without bias for all three files. There is some bias seen in the beginning from bases 1-10, though this could be due to adaptors.
- **Per sequence GC content**: The per sequence GC content for all three files demonstrated that the sequences have a normal distrubution of overall GC content signifying that the samples are likely not contaminated.
- **Per base N content**: This paramater shows that there was a low call of Ns (or unknown bases) in the sequences for all three files.
- **Sequence Length Distribution**: This parameter was flagged for all three files with an "!" meaning that it is not ideal but not detrimental. This metric looks at the size of the sequence fragments. Some sequencers keep this relatively uniform with a normal distribution, whereas others can output reads of varying lengths showing increasing distribution on the graphs.
- **Sequence Duplication Levels**: The sequence duplication levels were flagged with "!" for two of the files "AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq" and "AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq." This parameter counts the degree of duplication for every sequence in the library and creates a plot showing the relative number of sequences with different degrees of duplication. For these files the percent duplication shows that about 15% of the reads were duplicated around 10 times. This is okay, it is a low level of duplication that may indicate a very high level of coverage of the target sequence. If this number was high (above 20% for a large number of reads) it could signify enrichment bias (like PCR duplication).
- **Overrepresented sequences**: There were no overrepresented sequences for any of the files.
- **Adapter Content**: The adapter content for all three files indicated low to no adaptor content. 

## Trim reads with FastP 
Fastp will automatically trim adaptor content, which is not really necessary since Pumpkin's files all show low to no adaptor content from the fastqc. Fastp will also trim out reads with high N content. We will also use flags to trim the first 10 base pairs off the reads because these showed some bias in the per base sequence content (however low), and we will filter anything with an average quality score of 28 out and anything that is less than 20bp in length.
- `fastp -i AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq -o R163.fq -f 10 --qualified_quality_phred 28 --length_required 20`

No adapter detected

Read1 before filtering:
total reads: 15850172,
total bases: 2388224481,
Q20 bases: 2317822638(97.0521%),
Q30 bases: 2205087463(92.3317%)

Read1 after filtering:
total reads: 15392560,
total bases: 2165270001,
Q20 bases: 2119478868(97.8852%),
Q30 bases: 2023136059(93.4357%)

Filtering result:
reads passed filter: 15392560,
reads failed due to low quality: 455513,
reads failed due to too many N: 2099,
reads failed due to too short: 0,
reads with adapter trimmed: 0,
bases trimmed due to adapters: 0

- `fastp -i AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq -o R160.fq -f 10 --qualified_quality_phred 28 --length_required 20`

No adapter detected

Read1 before filtering:
total reads: 430051432,
total bases: 64402029541,
Q20 bases: 62491745953(97.0338%),
Q30 bases: 59427089941(92.2752%)

Read1 after filtering:
total reads: 417529051,
total bases: 58342047891,
Q20 bases: 57096870097(97.8657%),
Q30 bases: 54482802212(93.3851%)

Filtering result:
reads passed filter: 417529051
reads failed due to low quality: 12471474
reads failed due to too many N: 50907
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

- `fastp -i AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq -o R170.fq -f 10 --qualified_quality_phred 28 --length_required 20`

No adapter detected

Read1 before filtering:
total reads: 415960080,
total bases: 62751336061,
Q20 bases: 60405119840(96.2611%),
Q30 bases: 56854227696(90.6024%)

Read1 after filtering:
total reads: 400922995,
total bases: 56472690729,
Q20 bases: 54913988356(97.2399%),
Q30 bases: 51889157727(91.8836%)

Filtering result:
reads passed filter: 400922995,
reads failed due to low quality: 15019183,
reads failed due to too many N: 17902,
reads failed due to too short: 0,
reads with adapter trimmed: 0,
bases trimmed due to adapters: 0

## Check quality again with fastqc after trimming 
Run fastqc as follows and download the output .html files with sftp.
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Documents\Pumpkin`
- `fastqc R163.fq` `get R163_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.
- `fastqc R160.fq` `get R160_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.
- `fastqc R170.fq` `get R170_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.

## Align Pumpkin’s DNA to the domestic cat reference genome Felis catus 9.0 using the Burrows Wheeler Aligner and the default parameters 
Download the reference genome.
Go to: https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5 and download the reference fasta assembly. Upload this to the server with sftp.
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Downloads`
- `put felis_catus9.0.tar`
- `tar -xvf felis_catus9.0.tar`
- `cd ncbi-genomes-2022-04-21`
- `gunzip GCA_000181335.5_Felis_catus_9.0_genomic.fna`

Then make an index of the reference genome
- `bwa index -a bwtsw GCA_000181335.5_Felis_catus_9.0_genomic.fna FelisCatus`

Align the reads
- `bwa mem GCA_000181335.5_Felis_catus_9.0_genomic.fna ../R163.fq > aln-R163.sam`
- `bwa mem GCA_000181335.5_Felis_catus_9.0_genomic.fna ../R160.fq > aln-R160.sam`
- `bwa mem GCA_000181335.5_Felis_catus_9.0_genomic.fna ../R170.fq > aln-R170.sam`

This took a really long time, I'm not sure if the commands actually finished running for R160 or R170, I needed to go off WiFi while running both and so the terminal got halted. I will go with what was aligned just for the sake of time! --I could increase the threads to make these commands run faster but since the files are very large, I will not bother with it --we still have one file to work with!

## Manipulating the files with Samtools and removing duplicate reads
For further analysis, the SAM files from the alignments must be converted to BAM files. 
- `samtools view -S -b aln-R163.sam > aln-R163.bam`
- `samtools view -S -b aln-R160.sam > aln-R160.bam` --> Parse error at line ___  error reading file "aln-R160.sam": No such file or directory --must be due to the alignment not completing
- `samtools view -S -b aln-R170.sam > aln-R170.bam` --> Parse error at line ___ error reading file "aln-R170.sam": No such file or directory --must be due to the alignment not completing

Since the parse errors have occured in the files that did not finish aligning to the reference genome, we will just continue analyzing the one file that did finish (R163).

The alignments produced are in random order with respect to their position in the reference genome. We want to be able to call for variants so we must manipulate these BAM files so that the alignments occur in an order positionally based upon their alignment coordinates on each chromosome.
- `samtools sort aln-R163.bam -o aln-R163.sorted.bam`

Now we can index the genome sorted BAM files to allows us to extract alignments overlapping particular genomic regions. This allows us to quickly display alignments in each genomic region and is required by some genome viewers.
- `samtools index aln-R163.sorted.bam`

When we checked the quality of the reads with Fastqc we saw that the sequence duplication levels were low for the files (and were not even flagged for file R163), so we will not filter for any duplicates (which could be done with Samblaster). 

## SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20
Next I will use FreeBays, a Bayesian genetic variant detector, to find small polymorphisms (such as SNPs and indels) in Pumpkin's DNA. This program will use the cat reference genome and the sorted bam file and output a variant call format (vcf) file.

- `freebayes -f GCA_000181335.5_Felis_catus_9.0_genomic.fna aln-R163.sorted.bam > aln-R163var.vcf`

## SnpEff can then be used to detect variant effects. 
The vcf or variant call format files list all the differences between Pumpkin's DNA and the reference genome. To learn more about these variants other than their genetic coordinates (such as if they are in a gene, an exon, cause a protein coding change) we need to annotate the vcf file. SnpEff is a program that can do this for us.

- `conda install -c bioconda snpeff`
- `snpEff Felis_catus_9.0.99 aln-R163var.vcf > EFF.R163.vcf`

The vcf file can be looked at by using `cat EFF.R163.vcf`. When this is typed, the vcf file is printed to the terminal screen and the constant printing can be interrputed with `ctrl c`. The columns in this vcf file show that the annotations do not seem to be right, as the chromosome numbers are not correct. Cats have 38 chromosomes A1-A3, B1-B4, C1-C2, D1-D4, E1-E3, F1-F2, X, Y. None of these chromosomes were properly indentified. This error could be due to not knowing the nature of the fasta files. R163 was the only file that was small enough to be aligned to the cat reference genome. It could be Pumpkin's DNA, or it could be DNA used for his teeth microbiome analysis, or food analysis. 

Because of this, the following steps will use the vcf file provided by Basepaws that has proper alignment of Pumpkin's own cat DNA to the cat reference genome and proper annotations. 

## Filter the vcf data
Since the gvcf file provided by Basepaws is so large, it will not load into R and we can't make a manhattan plot of all the SNPs on all of Pumpkin's chromosomes. We need to filter the gvcf for what we need. I will filter for the cat chromosome D1 (10? Basepaws also doesn't have the cat chromosome numbers as A1-A3, B1-B4, C1-C2, D1-D4, E1-E3, F1-F2, X, Y...) because Basepaws indicated that Pumpkin had a SNP on this chromosome related to his health!

- `vcftools --vcf pumpkin_102.hard-filtered.gvcf --chr 10 --out filtVCF --recode`

## Visualize VCF data in R 
Next to visualize Pumpkin's SNPs we need to move to R Studio. this was really a nightmare because R kept freezing and I had to restart it each time it froze. When it did this, I had to reload the variables which took a long time. 

```
#load libraries
library(vcfR) 
library(tidyverse)

#load vcf file
BPvcf <- read.vcfR("~/FinalProject/filtVCF.recode.vcf", verbose = FALSE)
BPvcfT <- read.table("~/FinalProject/filtVCF.recode.vcf")
#load counts file
#vcf_counts <-  read.table("~/FinalProject/countsVCF.frq.count")

#get information from the vcfR file
chrom <- create.chromR(name="Chrom10", vcf=BPvcf)
#chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)
head(chrom)
pchrom <- plot(chrom)
qcchrom <- chromoqc(chrom, dp.alpha = 22)

#these graphs have shown us the quality of the vcf file, I probably could have filtered it more with vcftools but I wans't really sure what parameters to include
```
![pchrom](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/mgregoire/pchrom.png)

![qcchrom](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/mgregoire/qccrhom.png)

## Analyze data further in python
I am not a fan of R, so I decided to further analyze the data in a language I am more comfortable with --python. Below is the Python script I ran to analyze the vcf file further. I will also link to the Jupyter notebook [here](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/mgregoire/VCFanalysis.ipynb).
I have Anaconda Navigator installed on my computer and I used the python terminal to download a vcf viewer to my computer to use in Jupyter
-`conda install -c conda-forge scikit-allel`
Then I moved to Jupyter:
```
#load modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import io
import os
# import scikit-allel
import allel

#load the vcf file as a df with allel
df = allel.vcf_to_dataframe(r'C:\Users\mjgregoire\Documents\Pumpkin\filtVCF.recode.vcf', fields='*', alt_number=2)
df.head(5)
df.tail(5)

#count the number of SNPs in the dataframe
df.is_snp.value_counts() 
#this means Pumpkin has no snps on this chromosome!!! So maybe I was wrong in thinking the cat chromosome D1 is chromosome 10... that's okay! 
#if only I knew more about cat genetics! 

#we can still use this vcf of chromosome 10 to look at things like hetrozygosity and variants across position on the chromosome

#read in the vcf as a python dict
callset = allel.read_vcf(r'C:\Users\mjgregoire\Documents\Pumpkin\filtVCF.recode.vcf', fields='*')
#check what the attributes of the vcf dict are
sorted(callset.keys()) 
#‘samples’ array contains sample identifiers extracted from the header line in the VCF file
#arrays with keys beginning ‘variants/’ come from the fixed fields in the VCF file
#arrays with keys beginning ‘calldata/’ come from the sample fields in the VCF file

#make a genotype array based on the vcf file
gt = allel.GenotypeArray(callset['calldata/GT'])
gt
#count the number of heterozygous calls per variant
gt.count_het(axis=1)
#count the number times each allele (0=reference, 1=first alternate, 2=second alternate, etc.) is observed for each variant
ac = gt.count_alleles()
ac

#check the number of filtered reads that support each of the reported alleles
callset = allel.read_vcf(r'C:\Users\mjgregoire\Documents\Pumpkin\filtVCF.recode.vcf', fields=['variants/DP', 'calldata/DP'])
sorted(callset.keys())
callset['variants/DP']
callset['calldata/DP']

#get the positions
pos = allel.SortedIndex(callset['variants/POS'])
pos
#define a function to plot the density of the variants based on position
def plot_windowed_variant_density(pos, window_size, title=None):
    # setup windows 
    bins = np.arange(0, pos.max(), window_size)
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
#look at the graph
plot_windowed_variant_density(pos, window_size=100000, title='Variant density')
#save the graph as a figure
plt.savefig(r'C:\Users\mjgregoire\Documents\Pumpkin\variantDensity.png')

```
![Variant Density](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/mgregoire/variantDensity.png)


## Compare with the conclusions sent in Pumpkin's report!
From my data analysis, I learned that I don't know that much about cat genetics and that the process of analzying data from mammalian whole genome sequencing requires a LOT of storage and computer processing. I was able to make a workflow to process raw sequencing files from the fasta files to make filtered variant call format files that can be further analyzed in R or python. I wish that I could spend more time learning about what the files I obtained mean and parse through a bunch of different analyses. But I realize people spend lifetimes studying these fields and I've only spent a few weeks on this project! I can't wait to show Pumpkin what I was able to do with his DNA, but below let me just explain how detailed the analysis from Basepaws was.

Pumpkin's report from Basepaws was broken down into three main parts: a breed profile, a dental report, and his genetic report.

1. His breed profile indicated that he is mostly a western breed cat (68.82%). A large portion of his genome (17.7%) was found to be similar to Maine Coon cats (which explains why he's SO big!), Norwegian forest cats (7.88%), he also shares some DNA with Abyssinian cats (5.01%), persian cats (9.59%), and British shorthair (8%) among many others. Based on this Basepaws assumes that his blood type is likely type A based on his European and American descent.

2. The report also included a dental report which utilized DNA amplified from Pumpkin's mouth microbiome. They scored Pumpkin's risk of periodontal disease as low, his risk of tooth resporption as low, and his risk of bad breath as medium (which I actually can attest to being stinky!). This was assessed by looking at the sequences of bacterial DNA obtained from Pumpkin's spit on the swab. Out of 182 microbes that were detected, 44 are associated with bad breath! (Diaphorobacter sp., Acidovorax sp., and Pseudomonas denitirificans are Pumpkin's top three high risk microbes for bad breath.) They also characterized what Pumpkin is eating based on other eukaryotic DNA they found, and they said they found duck DNA. Pumpkin is mainly a seafood guy and he's only an indoor cat, so we think this is probably due to a filler in his dry food!

3. From his genetic report, he was found to be a carrier of one trait marker (dilute coat color, del(T) in gene MLPH), he was cleared of 31 genetic markers for disease, but was found to be at risk for one genetic marker involved in hypertrophic cardiomyopathy (HCM), where in the gene MYBPC3 he has one copy of a C>G mutation. This gene is found on the cat chromosome D1 and that's why I chose that chromosome to analyze from the very large gvcf file (I guess I was wrong in thinking chromosome D1 was 10. I assumed this because the Basepaws vcf file only had chromosomes in numbers not like A1-D1, so I had counted from A1 to D1 to get the 10). HCM is the most common feline heart disease and is characterized by thickening of the heart's muscular walls and tachycardia (fast heart beat). The severity of the disease varies and if it is diagnosed early treatment is easier. It is good to know this and inform Pumpkin's vet that he is a carrier for this condition! 

To recap, I never knew that the analysis behind all of Pumpkin's DNA would take so much work. I spent a couple weeks analyzing his DNA for this project and did not get anywhere near the extent of analysis that Basepaws did. It is truly amazing what can be done when you know what to look for and have the right equipment and an established pipeline. Kudos to all the wonderful people at Basepaws who made both analyzing Pumpkin's DNA for this project possible, and also knowing all about his health :)

## References
1. Basepaws: https://basepaws.com/
2. Fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
3. FastP: https://github.com/OpenGene/fastp
4. BWA: http://bio-bwa.sourceforge.net/
5. Samtools: http://www.htslib.org/
6. FreeBayes: https://github.com/freebayes/freebayes
7. SnpEff: https://pcingola.github.io/SnpEff/
8. VCF tools: http://vcftools.sourceforge.net/
9. vcfR: https://cran.r-project.org/web/packages/vcfR/index.html
10. Scikit-allel: https://scikit-allel.readthedocs.io/en/stable/
