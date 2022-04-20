# BIO594 Final Project Outline
Michelle Gregoire

## Option 1: Pumpkin’s DNA
This past year my family got our cat’s (Pumpkin’s) DNA sequenced! Whole genome sequencing was performed by Basepaws, and we recently received a flash drive with the raw sequence (fastq) files. For my final project I would like to analyze Pumpkin’s DNA if it is possible!

From the fastq files I would perform a fastqc to check for quality. The reads would then be trimmed with Trimmomatic or Cutadapt for anything below a quality score of 28, any reads with less than 20 base pairs, or any adaptor sequences. The quality after trimming would be checked once again with fastqc. Pumpkin’s DNA would then be mapped to the domestic cat reference genome Felis catus 9.0 (https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5) using the Burrows Wheeler Aligner (BWA if the reads are <100bp long, or BWA-MEM if the reads are 70bp to 1Mbp) and the default parameters. After alignment, Samblaster can then be used to identify and remove any duplicate reads. SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20. SnpEff can then be used to detect variant effects. 

## References:
- Buckley RM, Davis BW, Brashear WA, Farias FHG, Kuroki K, Graves T, et al. (2020) A new domestic cat genome assembly based on long sequence reads empowers feline genomic medicine and identifies a novel gene for dwarfism. PLoS Genet 16(10): e1008926. https://doi.org/10.1371/journal.pgen.1008926 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008926
- Gandolfi B, Grahn RA, Creighton EK, et al. COLQ variant associated with Devon Rex and Sphynx feline hereditary myopathy. Anim Genet. 2015;46(6):711-715. doi:10.1111/age.12350  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4637250/
- Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009;25(14):1754-1760. doi:10.1093/bioinformatics/btp324 ttps://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/
- Lyons LA, Creighton EK, Alhaddad H, et al. Whole genome sequencing in cats, identifies new models for blindness in AIPL1 and somite segmentation in HES7. BMC Genomics. 2016;17:265. Published 2016 Mar 31. doi:10.1186/s12864-016-2595-4 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4815086/)
- Xu X, Sun X, Hu XS, et al. Whole Genome Sequencing Identifies a Missense Mutation in HES7 Associated with Short Tails in Asian Domestic Cats. Sci Rep. 2016;6:31583. Published 2016 Aug 25. doi:10.1038/srep31583 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4997960/)

