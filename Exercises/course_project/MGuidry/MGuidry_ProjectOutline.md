# BIO 594 Bioinformatics Final Project Outline
**Megan Guidry**

# Project Title: Weighted Gene Co-Expression Network Analysis on Larval Oyster RNASeq Data
For the final project, I will do a weighted gene co-expression network analysis (WGCNA) on an existing, processed RNASeq data set.

## Brief Study Background
This study is focused on how independent and simultaneous stressors impact the function of early life stage oysters. The experiment exposed pools of *Crassostrea virginica* oyster larvae (trocophore stage) to treatment conditions for 24 hours. The treatments included: Control (ambient seawater), Coastal Acidification (2800 µatm pCO2 seawater), Sewage Effluent (5% treated sewage water in seawater), and a combined Coastal Acidification & Sewage Effluent (2800 µatm pCO2 & 5% sewage water seawater). These treatments are abreviated CO, CA, SE, and CASE, respectively. Larval pool subsamples were taken and flash frozen before and immediately following exposure to preserve DNA and RNA. 

The data for this project originate from 1 replicate block of CON, CA, SE and CASE treatments. There are 14 samples (3 CON, 4 CA, 3 SE, 4 CASE). Extracted RNA from this replicate block was with the KAPA Biosystems Stranded mRNA Seq kit and sequenced on one Full HiSeq lane from Novogene. 

## Status of project
The RNASeq Analysis (QC, filtering, aligning reads to oyster genome, assembling alignments to transcripts, counting transcripts, DESeq2 pipeline, and gene ontology analysis) were completed by Maggie Schedl. That analysis and more detailed project background can be found in [Maggie's Larval-Oyster-CASE-RNA Github repository](https://github.com/mguid73/Larval-Oyster-CASE-RNA). 🙌

## Applying a [Weighted Gene Co-Expression Network Analysis (WGCNA)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)
This analysis identifies clusters (called modules) of highly correlated genes, describes module characteristics (eigengene values and hub genes), relates modules to each other and sample metadata, and measures strength of module membership. Network correlations like this can be powerful tools to identify potential genes of interest after exposure to stress. 

