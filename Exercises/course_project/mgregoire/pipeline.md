# Analyzing Pumpkin's DNA 
## From the fastq files perform a fastqc to check for quality
Make sure that fastqc is is installed on your server. I am using the URI HPC Andromeda
- mkdir BIO594 | cd BIO594 | mkdir bioconda
- curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
- sh Miniconda3-latest-Linux-x86_64.sh
- conda config --add channels defaults
- conda config --add channels bioconda
- conda config --add channels conda-forge

Unzip the files
- gunzip AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq.gz
- gunzip AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq.gz
- gunzip AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq.gz

Run fastqc
- fastqc AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq
- fastqc AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq
- fastqc AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq

## Trim reads with Trimmomatic or Cutadapt for anything below a quality score of 28, any reads with less than 20 base pairs, or any adaptor sequences
- trimmomatic SE R163.fq AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20
- trimmomatic SE R160.fq AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20
- trimmomatic SE R170.fq AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20

## Check quality after trimming with fastqc
- fastqc R163.fq 
- fastqc R160.fq
- fastqc R170.fq

## Align Pumpkin’s DNA to the domestic cat reference genome Felis catus 9.0 using the Burrows Wheeler Aligner and the default parameters 
Download the reference genome
- wget https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5
First make an index of the reference genome
- bwa index [-a bwtsw|is] GCA_000181335.5.fasta FelisCatus
Align the reads
-

## After alignment, Samblaster can then be used to identify and remove any duplicate reads
- samtools view -h samp.bam | samblaster --ignoreUnmated [-e] --maxReadLength 100000 [-s samp.split.sam] [-u samp.umc.fasta] | samtools view -Sb - > samp.out.bam
- samtools view -h samp.bam | samblaster --ignoreUnmated -a [-e] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null

## SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20
- freebayes -f ref.fa aln.bam >var.vcf
- env/bin/platypus callVariants --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf

## SnpEff can then be used to detect variant effects. 
