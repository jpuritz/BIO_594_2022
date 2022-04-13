# Project Title: 
# Caribbean Scale Population Structure of Deep & Shallow Lineages of the Octocoral *Eunicea flexuosa*

## Project Rational
The Caribbean wide octocoral *Eunicea flexuosa* has been shown to posses morphologically and genetically distinct sympatric lineages that appear to segregate by depth. Work thus far has shown that there is minimal mixing of these lineages within islands. The proportion of individuals originating from deep or shallow lineages changes predictably with depth, and evidence suggests that if a hybrid zone exists it is small. The prevailing hypothesis is that these two lineages likely represent recently diverged cryptic species where the exchange of genetic material is minimal. The genetic barrier is thought to be maintained through immigrant inviability, and unsuccessful hybrids. Further evidence suggests that the genetic barrier is reinforced by a slight (hour scale) difference in gamete release during spawning. 

It is not known how these populations relate to one another across islands. This project seeks to explore the population structure of these two cryptic lineages across the Caribbean basin. If deep and shallow lineages are more like their respective lineage across islands compared to their opposing sympatric lineage that would be strong evidence in support of two cryptic species of *Eunicea flexuosa*.

## Project Plan & Goals
1. Creation of SAM/BAM files from demultiplexed ddRAD samples using `BWA` and `SAMtools`

2. Creation of *de novo* reference using `dDocent` for comparison of individuals across and within islands from each lineage 

3.  Quality filtering of SNP reads using `dDocent` and `VCF tools`

------------------------------------------------------------------
4.  Detection of outlier SNPs using `Bayscan`and removal of outlier SNPs through subsetting

5.  Creation of Fst population structure statistics from ddRAD samples of *Eunicea flexuosa* taken across the Caribbean basin (could test with `hierfstat` in `R`)

6. Create an Fst structure plot of lineage by locality (island) (`OutFlank` may be a good place to start)

7. Generation of PCA figures to visualize population structure (could use `PCA adapt` or `adegenet` &
`vcfR`)
------------------------------------------------------------------

## Project Details

### Data

Individual octocorals were morphologically identified to lineage and sampled across various depths (<4m – ~30m) throughout the Caribbean basin (Puerto Rico, Florida, Bahamas, Panama, Honduras, Guadalup, Dominican Republic, & Curacao). Libraries were prepped using a ddRAD protocol where a combination of two unique adapters identifies individuals.

Data files are demultiplexed gziped fastq files of the forward and reverse reads of each individual. Examples `bc10_pa_1-AGTCATCGG_f_pa.fq.gz`, `bc10_pa_1-AGTCATCGG_r_pa.fq.gz` 

### Project Environments
* Objectives 1–3, and parts of 4 will be performed on the URI Bluewaves server where the data is stored.
* Objectives 5–6, and parts of 4 will be performed in R