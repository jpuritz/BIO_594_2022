# BIO 594 Course Project Outline

# Hybrid assembly of cDNA and EecSeq reads using larval oyster samples under CASE stress

## Problem Statement

EecSeq is a new tool for evolutionary biologists to explore questions concerning genes under selection. By reducing the genomic information captured we can greatly expand the number of samples a sequencing study design can accommodate. Additionally EecSeq has no need for a reference to design probes opening up the potential to study non-model or obscure model organisms, but without a reference EecSeq reads must be assembled through a de novo bioinformatic framework. The EecSeq protocol creates cDNA libraries of high quality that could be sequenced and guide EecSeq assembly through a de novo transcriptome using these previously unsequenced cDNA libraries. This project will explore de novo assembly of RNAseq reads into an annotated de novo transcriptome that will act as a sudo-reference for EecSeq read assembly. 

## Methods

### Sources of Data

* RNAseq reads from CASE study
* EecSeq reads from CASE study
* Experimental design map
    * including library preparation
    * exposure condition

### Bioinformatic approach

![](BIO594_project_assembly_outline.png)

### Read Filtering

Both cDNA and gDNA Reads will be symbolically linked to respective working directories within the home directory of the user on KITT. [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) will assess each read pool and [Multiqc](https://multiqc.info/) will be used to visualize difference in read quality and content, specifically for per sequence quality scores, N content, and overrepresetned sequences. Reads will then be processed using fastp to trim adapters, low quality bases, and problematic reads. [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Multiqc](https://multiqc.info/) will be used again to assess read quality and visualize.

There are three approaches to consider for read normalization which differ in efficacy between cDNA and gDNA reads. cDNA reads benefit from in silico read normalization (typically through [Trinity's](https://github.com/trinityrnaseq/trinityrnaseq) insilico_read_normalization.pl script) where gDNA reads are filtered through kmer identifcation or unique read filtering. In our case we will apply [Trinity's](https://github.com/trinityrnaseq/trinityrnaseq) read normalization at 50 x coverage to the cDNA reads and [dDocent](https://www.ddocent.com/) unique read filtering to our gDNA reads.

### De novo assembly

For de novo assembly of both read pools we will be applying three separate de novo assembly programs: [Trinity](https://github.com/trinityrnaseq/trinityrnaseq), [Trans-ABySS](https://github.com/bcgsc/transabyss), and [Oases](https://github.com/dzerbino/oases). One important factor to consider in these approaches is alteration to kmer length assembly. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) limits kmer length to <= 32 base pairs, while [Trans-ABySS](https://github.com/bcgsc/transabyss), and [Oases](https://github.com/dzerbino/oases) can be set higher but no greater than 2/3 read length. Any kmer options set higher than 2/3 read length can result in fewer contigs. Kmer lengths that are too short can lead to many short contigs, which many are partial gene representations. For both [Trans-ABySS](https://github.com/bcgsc/transabyss), and [Oases](https://github.com/dzerbino/oases) we will iterate kmer lengths starting at 25 moving up by 12 till 2/3 read length to determine the optimal de novo assembly based on assembly statistics.

### Clustering

[CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) will be used to cluster the initital de novo assemblies into a draft assembly. Clustering of de novo assemblies reduces redundancy and inproves downstream analysis by removing multi-copy contigs. More specifics on thresholds and options to come.

### Mapping

[BWA](http://bio-bwa.sourceforge.net/+) will be our go to mapper through this project. There are two different parts to mapping in this project. The first is mapping filtered gDNA reads to the transcriptome in order to guide Exome assembly. The second is to remap gDNA reads to the Exome assembly guided by the transcriptome assembly. [BWA](http://bio-bwa.sourceforge.net/+) creates an index from the reference first and them will use the PE-mapping BWA-MEM algorithm which will output a SAM file which can be converted on the fly into a lower footprint BAM file which is sorted. 

* [Salmon](https://github.com/COMBINE-lab/salmon)

### Assembly assessment

* [Trinity stats](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* [DETONATE](https://github.com/deweylab/detonate)
* [Transrate](https://github.com/blahah/transrate)
* [RNAquast](https://github.com/ablab/rnaquast)
* [BUSCO](https://busco.ezlab.org/)
* [SOURMASH?](https://sourmash.readthedocs.io/en/latest/)

### Variant calling analysis

* [VCFtools](https://vcftools.github.io/index.html)
* [PCA](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/EecSeq_Cvirginica_OutlierDetection.md)
* [DAPC](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/EecSeq_Cvirginica_OutlierDetection.md)
* [RDA](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/PopGen_SeaGen_Analyses/RedundancyAnalysis/RDA_Outlier_Hap.Rmd)

### References

