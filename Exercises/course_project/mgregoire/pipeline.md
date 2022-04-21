# Analyzing Pumpkin's DNA
## WGS Data Analysis Pipeline
![Pumpkin](https://github.com/jpuritz/BIO_594_2022/tree/main/Exercises/course_project/mgregoire/pumpkin_resize.jpg)
## Set up directories, softwares, and files 
Make sure that bioconda is is installed. I am using the Kitt server and bioconda has already been installed. Bioconda makes it easy to install many of the various programs that will be used in this pipeline (eg: fastqc).
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
- **Per base sequence quality**: The quality scores across all bases look good for all 3 fastq files! The reads start to drop off towards the tail ends (this is expected as the reads get longer), but the scores stay in the green zone the entire time hovering around 36-34 and only with error bars dropping to about 25 from 100-150bp. This means that the average of the reads had >99.9% accuracy.
- **Per sequence quality scores**: The quality score distribution across all sequences peak at 36 for all three files. This means that Q36 had more read numbers than other quality scores.
- **Per base sequence content**: This metric shows that the percentage of each of the four nucleotides (A, T, C, G) at each position across all reads showed normal random distributions without bias for all three files. There is some bias seen in the beginning from bases 1-10, though this could be due to adaptors.
- **Per sequence GC content**: The per sequence GC content for all three files demonstrated that the sequences have a normal distrubution of overall GC content signifying that the samples are likely not contaminated.
- **Per base N content**: This paramater shows that there was a low call of Ns (or unknown bases) in the sequences for all three files.
- **Sequence Length Distribution**: This parameter was flagged for all three files with an "!" meaning that it is not ideal but not detrimental. This metric looks at the size of the sequence fragments. Some sequencers keep this relatively uniform with a normal distribution, whereas others can output reads of varying lengths showing increasing distribution on the graphs.
- **Sequence Duplication Levels**: The sequence duplication levels were flagged with "!" for two of the files "AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq" and "AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq." This parameter counts the degree of duplication for every sequence in the library and creates a plot showing the relative number of sequences with different degrees of duplication. For these files the percent duplication shows that about 15% of the reads were duplicated around 10 times. This is okay, it is a low level of duplication that may indicate a very high level of coverage of the target sequence. If this number was high (above 20% for a large number of reads) it could signify enrichment bias (like PCR duplication).
- **Overrepresented sequences**: There were no overrepresented sequences for any of the files.
- **Adapter Content**: The adapter content for all three files indicated low to no adaptor content. 

## Trim reads with FastP 
Fastp will automatically trim adaptor content, which is not really necessary since Pumpkin's files all show low to no adaptor content from the fastqc. Fastp will also trim out reads with high N content. We will also use flags to trim the first 10 base pairs off the reads because these showed some bias in the per base sequence content (however low), and we will filter anything with an average quality score of 28 out and anything that is less than 20bp in length.
- `fastp -i AA.WQ.58.31210251600500.LP.713.F3.L2.R163.WGS.fastq -o R163.fq -f 10 --qualified_quality_phred 28 --length_required 20`

No adapter detected

Read1 before filtering:
total reads: 15850172,
total bases: 2388224481,
Q20 bases: 2317822638(97.0521%),
Q30 bases: 2205087463(92.3317%)

Read1 after filtering:
total reads: 15392560,
total bases: 2165270001,
Q20 bases: 2119478868(97.8852%),
Q30 bases: 2023136059(93.4357%)

Filtering result:
reads passed filter: 15392560,
reads failed due to low quality: 455513,
reads failed due to too many N: 2099,
reads failed due to too short: 0,
reads with adapter trimmed: 0,
bases trimmed due to adapters: 0

- `fastp -i AA.WQ.58.31210251600500.SP.299.B2.L2.R160.WGS.fastq -o R160.fq -f 10 --qualified_quality_phred 28 --length_required 20`

No adapter detected

Read1 before filtering:
total reads: 430051432,
total bases: 64402029541,
Q20 bases: 62491745953(97.0338%),
Q30 bases: 59427089941(92.2752%)

Read1 after filtering:
total reads: 417529051,
total bases: 58342047891,
Q20 bases: 57096870097(97.8657%),
Q30 bases: 54482802212(93.3851%)

Filtering result:
reads passed filter: 417529051
reads failed due to low quality: 12471474
reads failed due to too many N: 50907
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

- `fastp -i AA.WQ.58.31210251600500.SP.307.D1.L1.R170.WGS.fastq -o R170.fq -f 10 --qualified_quality_phred 28 --length_required 20`
No adapter detected

Read1 before filtering:
total reads: 415960080,
total bases: 62751336061,
Q20 bases: 60405119840(96.2611%),
Q30 bases: 56854227696(90.6024%)

Read1 after filtering:
total reads: 400922995,
total bases: 56472690729,
Q20 bases: 54913988356(97.2399%),
Q30 bases: 51889157727(91.8836%)

Filtering result:
reads passed filter: 400922995,
reads failed due to low quality: 15019183,
reads failed due to too many N: 17902,
reads failed due to too short: 0,
reads with adapter trimmed: 0,
bases trimmed due to adapters: 0

## Check quality again with fastqc after trimming 
Run fastqc as follows and download the output .html files with sftp.
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Documents\Pumpkin`
- `fastqc R163.fq` `get R163_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.
- `fastqc R160.fq` `get R160_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.
- `fastqc R170.fq` `get R170_fastqc.html`--> The trimming cleaned up the per base sequence quality and per base sequence content.

## Align Pumpkin’s DNA to the domestic cat reference genome Felis catus 9.0 using the Burrows Wheeler Aligner and the default parameters 
Download the reference genome
Go to: https://www.ncbi.nlm.nih.gov/assembly/GCA_000181335.5 and download the reference fasta assembly. Upload this to the server with sftp.
- `sftp -P {port number here} {username}@kitt.uri.edu`
- `cd FinalProject`
- `lcd Downloads`
- `put felis_catus9.0.tar`
- `tar -xvf felis_catus9.0.tar`
- `cd ncbi-genomes-2022-04-21`
- `gunzip GCA_000181335.5_Felis_catus_9.0_genomic.fna`
Then make an index of the reference genome
- `bwa index GCA_000181335.5_Felis_catus_9.0_genomic.fna FelisCatus`
Align the reads
- `bwa mem ref.fa ../R163.fq > aln-R163.sam`
- `bwa mem ref.fa ../R160.fq > aln-R160.sam`
- `bwa mem ref.fa ../R170.fq > aln-R170.sam`

## After alignment, Samblaster can then be used to identify and remove any duplicate reads
- `samtools view -h samp.bam | samblaster --ignoreUnmated [-e] --maxReadLength 100000 [-s samp.split.sam] [-u samp.umc.fasta] | samtools view -Sb - > samp.out.bam`
- `samtools view -h samp.bam | samblaster --ignoreUnmated -a [-e] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null`

## SNPs and small indels can be investigated with Samtools 1.2, Platypus, or FreeBayes, filtering out anything with a low quality of less than 20
- `freebayes -f ref.fa aln.bam >var.vcf`
- `env/bin/platypus callVariants --bamFiles=input.bam --refFile=reference.fa --output=variant_calls.vcf`

## SnpEff can then be used to detect variant effects. 
