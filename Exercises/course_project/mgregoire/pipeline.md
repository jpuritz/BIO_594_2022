# Analyzing Pumpkin's DNA 
![Pumpkin](https://github.com/jpuritz/BIO_594_2022/tree/main/Exercises/course_project/mgregoire/pumpkin_resize.jpg)
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
Run fastqc and download the .html output files
- `fastqc *.fastq`
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Documents\Pumpkin`
-  `get AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS_fastqc.html` | get AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS_fastqc.html | get AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS_fastqc.html`

The fastqc files from Pumpkin's DNA reveal the following:
- Per base sequence quality: The quality scores across all bases look good for all 3 fastq files! The reads start to drop off towards the tail ends (this is expected as the reads get longer), but the scores stay in the green zone the entire time hovering around 36-34 and only with error bars dropping to about 25 from 100-150bp.
- Per sequence quality scores: The quality score distribution across all sequences peak at 36 for all three files. This means that Q36 had more read numbers than other quality scores.
- Per base sequence content: This metric shows that the percentage of each of the four nucleotides (A, T, C, G) at each position across all reads showed normal random distributions without bias for all three files. There is some bias seen in the beginning from bases 1-10, though this could be due to adaptors.
- Per sequence GC content: The per sequence GC content for all three files demonstrated that the sequences have a normal distrubution of overall GC content signifying that the samples are likely not contaminated.
- Per base N content: This paramater shows that there was a low call of Ns (or unknown bases) in the sequences for all three files.
- Sequence Length Distribution: This parameter was flagged for all three files with an "!" meaning that it is not ideal but not detrimental. This metric looks at the size of the sequence fragments. Some sequencers keep this relatively uniform with a normal distribution, whereas others can output reads of varying lengths showing increasing distribution on the graphs.
- Sequence Duplication Levels: The sequence duplication levels were flagged with "!" for two of the files "AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq" and "AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq." This parameter counts the degree of duplication for every sequence in the library and creates a plot showing the relative number of sequences with different degrees of duplication. For these files the percent duplication shows that about 15% of the reads were duplicated around 10 times. This is okay, it is a low level of duplication that may indicate a very high level of coverage of the target sequence. If this number was high (above 20% for a large number of reads) it could signify enrichment bias (like PCR duplication).
- Overrepresented sequences: There were no overrepresented sequences for any of the files.
- Adapter Content: The adapter content for all three files indicated low to no adaptor content. 

## Trim reads with Trimmomatic for the first 10 base pairs, anything below an average quality score of 28, and any reads less than 20 base pairs
We will trim the first 10 base pairs off the reads because these showed some bias in the per base sequence content (however low). We will also filter anything with an average quality score of 28 out (though the average for all three files was around 36 so this is futile), and anything that is less than 20bp (though fastqc said the sequences were all between 31-151bp so this is also futile).
- `trimmomatic SE R163.fq AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq HEADCROP:10 AVGQUAL:28 MINLEN:20`
- `trimmomatic SE R160.fq AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq HEADCROP:10 AVGQUAL:28 MINLEN:20`
- `trimmomatic SE R170.fq AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq HEADCROP:10 AVGQUAL:28 MINLEN:20`

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
