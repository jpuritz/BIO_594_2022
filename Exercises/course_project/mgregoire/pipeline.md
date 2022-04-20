# From the fastq files perform a fastqc to check for quality

# Trim reads with Trimmomatic or Cutadapt for anything below a quality score of 28, any reads with less than 20 base pairs, or any adaptor sequences

# Check quality after trimming with fastqc

# Align Pumpkin’s DNA to the domestic cat reference genome Felis catus 9.0 (https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5) using the Burrows Wheeler Aligner (BWA if the reads are <100bp long, or BWA-MEM if the reads are 70bp to 1Mbp) and the default parameters. 

# After alignment, Samblaster can then be used to identify and remove any duplicate reads

# SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20

# SnpEff can then be used to detect variant effects. 
