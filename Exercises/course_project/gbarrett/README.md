# Genome-Wide Association NextFlow Pipeline 

Pipeline for obtaining Bayes Factor (BF), Population Differientiation (XtX), Fst z-scores, and latent factor mixed model (lfmm) z-scores. 

## Background

Genome-wide association and differientation methods are a valuable tool for honing an analysis to loci putatively under selection. To elicit a list of candidate loci, population and/or individual genotypes and environmental data need to be provided. In and effort to expedite this process I have glued together several bash scripts into nextflow. With the hope that through channel manipulation I can perform environmental GWAS tests in parallel. 

## Usage

`nextflow run main.nf`

## ToDo

    1. write docker container to run MINOTAUR R package
    2. add MINOTAUR
    2. High SNPs in Consensus vcf: VCF stats potentially due to bcftools norm 
    4. switch to lfmm2
    5. add RDA


## Final Comments

Due to time constraints I was unable to combine gwas statistics described above into a composite statistic for manhattan plot visualization. I also was unable to include RDA analysis and visualizations.