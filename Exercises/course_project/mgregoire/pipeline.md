# Analyzing Pumpkin's DNA 
![Pumpkin](https://github.com/jpuritz/BIO_594_2022/tree/main/Exercises/course_project/mgregoire/pumpkin_resize.jpg?raw=true)
## Set up your directories, softwares, and files 
Make sure that bioconda is is installed on your server. I am using the Kitt server and bioconda has already been installed. Bioconda makes it easy to install many of the various programs that will be used in this pipeline (eg: fastqc).
Make a folder for the project and cd into it. 
- `mkdir FinalProject | cd FinalProject`

Since the files were supplied to me on a flashdrive from Basepaws, I sftp'd the files onto the kitt server with the following commands:
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd D: | lls`
- `put AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq.gz | put AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq.gz | put AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq.gz | put pumpkin_102.hard-filtered.gvcf | put pumpkin_102.hard-filtered.gvcf.gz.tbi`

Unzip the files
- `gunzip AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq.gz`
- `gunzip AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq.gz`
- `gunzip AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq.gz`

## Run fastqc on the files to check for quality
Run fastqc
- `fastqc *.fastq`

## Trim reads with Trimmomatic or Cutadapt for any adaptor sequences, anything below an average quality score of 28, and any reads less than 20 base pairs
- `trimmomatic SE R163.fq AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20`
- `trimmomatic SE R160.fq AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20`
- `trimmomatic SE R170.fq AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq ILLUMINACLIP [[HEADCROP:5 (change this if need first 5 cut)]]] AVGQUAL 28 MINLEN:20`

## Check quality after trimming with fastqc
- `fastqc *.fq`

## Align Pumpkin’s DNA to the domestic cat reference genome Felis catus 9.0 using the Burrows Wheeler Aligner and the default parameters 
Download the reference genome
- `wget https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5`
First make an index of the reference genome
- `bwa index [-a bwtsw|is] GCA_000181335.5.fasta FelisCatus`
Align the reads
- `bwa mem ref.fa reads.fq > aln-se.sam`

## After alignment, Samblaster can then be used to identify and remove any duplicate reads
- `samtools view -h samp.bam | samblaster --ignoreUnmated [-e] --maxReadLength 100000 [-s samp.split.sam] [-u samp.umc.fasta] | samtools view -Sb - > samp.out.bam`
- `samtools view -h samp.bam | samblaster --ignoreUnmated -a [-e] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null`

## SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20
- `freebayes -f ref.fa aln.bam >var.vcf`
- `env/bin/platypus callVariants --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf`

## SnpEff can then be used to detect variant effects. 
