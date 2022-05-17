# Genome-Wide Association NextFlow Pipeline 

Pipeline for obtaining Bayes Factor (BF), Population Differientiation (XtX), Fst z-scores, and latent factor mixed model (lfmm) z-scores. 

## Background

Genome-wide association and differientation methods are a valuable tool for honing an analysis to loci putatively under selection. To elicit a list of candidate loci, population and/or individual genotypes and environmental data need to be provided. In and effort to expedite this process I have glued together several bash scripts into nextflow. With the hope that through channel manipulation I can perform environmental GWAS tests in parallel. 

## Usage

`nextflow run main.nf` 

## ToDo

    1. write docker container 
    2. High SNPs in Consensus vcf: VCF stats
    3. switch to lfmm2
    4. add RDA
    5. add MINOTAUR