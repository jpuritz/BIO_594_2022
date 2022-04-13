# Week 10 Eecseq, PCA, NMDs, and Pooled data

## Goals

[Objective from github](https://github.com/jpuritz/BIO_594_2022/tree/main/Exercises/Week10)

### EecSeq

1. To explore a small EecSeq data set
2. To learn how to use bedtools and samtools to explore coverage across genomic intervals and annotations
3. To plot the coverage of EecSeq pools across a gene of interest

## Setup Data and Environment

Loading ANGSD biconda environment

```bash
mkdir wk10
cd wk10

#Creating link to all in data folder | each file is a symbolic link
ln -s /home/BIO594/Exercises/Week_10/data/* .

#Making a file from a list of all .bam files

ls *.bam > all.files
```
>We have two file types in these links .bam and .bam.bai `all.files` has both

## Start of EecSeq Protocol

```bash
mkdir eecseq
cd eecseq
#creating environment
mamba create -n eecseq ddocent
source activate eecseq 

#Linking data from class directory

ln -s /home/BIO594/Exercises/Week_10/EecSeq/* .
```
> We have F and R.fq.gz files, reference fasta files, .gff3, and exon.cov.stats files

Overall there are three replicated capture pools EC 2, 4, & 7., a chromosome 1 file, and reference files.

The data using exome capture has already been trimmed and only contains reads from chromosome 1 of the eastern oyster genome.

Moving on to BWA and sorting with samtools.

[bwa manual ref pages](http://bio-bwa.sourceforge.net/bwa.shtml)

```bash
# EC2
bwa mem reference.fasta EC_2.F.fq.gz EC_2.R.fq.gz -t 8 -a -M -B 3 -O 5 -R "@RG\tID:EC_2\tSM:EC_2\tPL:Illumina" 2> bwa.EC_2.log | samtools view -@4 -q 1 -SbT reference.fasta - > EC_2.bam
# EC4
bwa mem reference.fasta EC_4.F.fq.gz EC_4.R.fq.gz -t 8 -a -M -B 3 -O 5 -R "@RG\tID:EC_4\tSM:EC_4\tPL:Illumina" 2> bwa.EC_4.log | samtools view -@4 -q 1 -SbT reference.fasta - > EC_4.bam
# EC 7
bwa mem reference.fasta EC_7.F.fq.gz EC_7.R.fq.gz -t 8 -a -M -B 3 -O 5 -R "@RG\tID:EC_7\tSM:EC_7\tPL:Illumina" 2> bwa.EC_7.log | samtools view -@4 -q 1 -SbT reference.fasta - > EC_7.bam
# EC 2
samtools sort -@8 EC_2.bam -o EC_2.bam && samtools index EC_2.bam 
# EC 4
samtools sort -@8 EC_4.bam -o EC_4.bam && samtools index EC_4.bam 
# EC 7
samtools sort -@8 EC_7.bam -o EC_7.bam && samtools index EC_7.bam 
```
Breaking down `bwa` 
* mem - specifies BWA-MEM method
* .fasta - reference file
* .fq.gz - input files
* -t - number of threads
* -a - output all found alignments for single-end or unpaired paired-end reads
* -M - mark shorter split hits as secondary
* -B - Mismatch penalty
* -O - Gap open penalty
* -R - Complete read group header line `/t`. Read group ID will be attached to every read in output An example is ’@RG\tID:foo\tSM:bar’

> Come back and finish the annotations below - 13-apr-2022

Breaking down sametools sort
* -@ - number of threads
* input file
* -o 
----------------------------------

### Marking duplicate reads

(Picard tools documentation)[https://broadinstitute.github.io/picard/]

```bash
# loading picard tools
wget https://github.com/broadinstitute/picard/releases/download/2.17.8/picard.jar

# Marking PCR duplicates EC 2

java -Xms4g -jar picard.jar MarkDuplicatesWithMateCigar I=EC_2.bam O=EC_2md.bam M=EC_2_dup_metrics.txt MINIMUM_DISTANCE=300
```
> Come back and annotate these flags 13-apr-2022

> Out put matches Jon's 250326 marking records, and 14208 optical duplicate clusters

```bash
# repeating for EC 4
java -Xms4g -jar picard.jar MarkDuplicatesWithMateCigar I=EC_4.bam O=EC_4md.bam M=EC_4_dup_metrics.txt MINIMUM_DISTANCE=300

# repeating for EC 7
java -Xms4g -jar picard.jar MarkDuplicatesWithMateCigar I=EC_7.bam O=EC_7md.bam M=EC_7_dup_metrics.txt MINIMUM_DISTANCE=300
```
> EC 4: Output does not match Jon's 215837 marking records (diff -1), 11009 optical duplicat clusters (diff -1).

> EC 7: Output does not match Jon's 196122 (diff -7) marking records, 9913  optical duplicate clusters (diff +1)

### Filtering out duplicates, and secondary alignments

Mappings with quality scores less than 10 will be removed and reads with more than 80bp will be clipped

```bash
# EC 2
samtools view -@8 -h -F 0x100 -q 10 -F 0x400 EC_2md.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@8 -b > EC_2.F.bam 
# EC 4
samtools view -@8 -h -F 0x100 -q 10 -F 0x400 EC_4md.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@8 -b > EC_4.F.bam 
# EC 7
samtools view -@8 -h -F 0x100 -q 10 -F 0x400 EC_7md.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@8 -b > EC_7.F.bam
```
> Go back and annotate flags - 13-apr-2022

> Once the semester ends find resources and time to understand `mawk` code

#### Verificaton of filtered reads

Commands utilize subshells allowing two instances of samptools to run "in parallel". We want to see two number such that the 2nd is less than the first

```bash
# "Paste the output of each shell (<), assuming -c means count.

# EC 2
paste <(samtools view -c EC_2md.bam) <(samtools view -c EC_2.F.bam )
# EC 4
paste <(samtools view -c EC_4md.bam) <(samtools view -c EC_4.F.bam )
# EC7
paste <(samtools view -c EC_7md.bam) <(samtools view -c EC_7.F.bam )
```
>EC2: Out put differs from Jon's but shows a reduction in reads 4452400 vs. 4046543

>EC4: Out put differs from Jon's but shows a reduction in reads 3640300 vs. 3287915

>EC7:Out put differs from Jon's but shows a reduction in reads 3642646 vs. 3316261

#### Exploratory analysis: depth per bp along refernce

samtools does this well

```bash
# EC 2
samtools depth -aa EC_2.F.bam > EC_2.genome.depth
# EC 4
samtools depth -aa EC_4.F.bam > EC_4.genome.depth
#EC 7
samtools depth -aa EC_7.F.bam > EC_7.genome.depth

# Checking out put

#EC 2
head EC_2.genome.depth
#EC 4
head EC_4.genome.depth
#EC 7
head EC_7.genome.depth

```

> Overall it looks like EC 4 has  the best coverage per bp per chromosome (at least for the first 10).

#### Subsetting to specific regions of chromosom 1 that influence heat shock proteins

```bash
#EC 2
mawk '$2 > 32736205 && $2 < 32866205' EC_2.genome.depth > EC_2.graph.depth
#EC 4 
mawk '$2 > 32736205 && $2 < 32866205' EC_4.genome.depth > EC_4.graph.depth
#EC 7 
mawk '$2 > 32736205 && $2 < 32866205' EC_7.genome.depth > EC_7.graph.depth 
```
# Will need to come back to this assignment later