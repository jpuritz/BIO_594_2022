Code for my bio594 project. Analysis done on Andromeda 

Proposal [here](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Ashey_BIO594_FinalProjectProposal.md)

idk who said this but wise words: "the computer is never wrong, it is you who is wrong"

## Table of Contents 

I. [Set-up](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#set-up-directories-and-data)

II. [Quality check of raw reads](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#quality-check)

III) [Trim raw reads](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#trim)

IV) [Quality check of trimmed reads](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#quality-check-trimmed-reads)

V) [Align reads](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#align-reads-using-different-read-aligners)
	- A) [Align against genome - STAR](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#a-align-against-genome---star)
		- i) [A.cervicornis genome index](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#i-acerv-genome-index)
		- ii) [P.acuta genome index](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#ii-pacuta-genome-index)
		- iii) [A.cervicornis alignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#iii-acerv-alignment)
		- iv) [Pacuta alignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#iv-pacuta-alignment)
	- B) [Align against de novo transcriptome - Trinity & Bowtie](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#b-align-against-transcriptome---trinity--bowtie2)
		- i) [A.cervicornis de novo transcriptome Trinity](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#i-acerv-trinity)
		- ii) [P.acuta de novo transcriptome Trinity](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#ii-pacuta-trinity)
		- iii) [A.cervicornis Bowtie alignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#iii-acerv-bowtie)
		- iv) [P.acuta Bowtie alignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#iv-pacuta-bowtie)
	- C) [Pseudoalignment - Kallisto](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#c-pseudoalignment---kallisto)
		- i) [A.cervicornis pseudoalignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#i-acerv-pseudoalignment)
		- ii) [P.acuta pseudoalignment](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#ii-pacuta-pseudoalignment)

VI) [Alignment comparisons](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#vi-compare)

VII) [References](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Workflow.md#vii-references)


### I) Set up directories and data 

Make directory for project

```
cd /data/putnamlab/jillashey/
mkdir BIO594_FinalProject
cd BIO594_FinalProject
```
Create folders inside 

```
mkdir data QC STAR Bowtie2 Kallisto Trinity 
cd data
mkdir raw trimmed
```

Make links to Acerv and Pacuta data 

```
cd raw

# Pacuta
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/1_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/2_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/4_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/11_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/28_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/31_2.fastq . 
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/35_2.fastq . 
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/36_2.fastq . 
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/38_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/39_2.fastq . 
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/41_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/42_2.fastq .
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/47_2.fastq .

# Acervicornis 
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/raw/*Ac* .
# I removed any Acerv data with an *_1.fastq, as I am just using the *_2.fastq files
```

Check number of files to make sure everything was linked properly

```
ls | wc -l
28
```

Count number of reads per file. Some reads have @HISEQ as header, some reads have @HWI

```
zgrep -c "@HISEQ" *.fastq > HISEQ_raw_length.txt

zgrep -c "HWI" *.fastq > HWI_raw_length.txt
```

### II) Quality check

Quality checking will be done using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info). 

Make folders for raw and trimmed QC information

```
mkdir raw trimmed
cd raw 
```

Run QC on raw reads 

```
nano fastqc_raw.sh 

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="fastqc_out_raw_error"
#SBATCH --output="fastqc_out_raw"

echo "START"; date

module load FastQC/0.11.9-Java-11 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

base="/data/putnamlab/jillashey/BIO594_FinalProject"

# Run fastqc on raw files 
for file in "$base"/data/raw/*.fastq
do
fastqc $file --outdir "$base"/QC/raw/
done

# Compile with MultiQC
cd "$base"/QC/raw/
multiqc *fastqc.zip -o .

echo "STOP"; date

sbatch fastqc_raw.sh 
```

Submitted batch job 129345

Copy MultiQC files onto computer and look at plots

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/raw_fastqc_per_base_sequence_quality_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/raw_fastqc_per_sequence_gc_content_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/raw_fastqc_adapter_content_plot.png)

I didn't add all the plots, just the ones that I thought were most important. Sequences were either 50 or 125 bp long. Quality scores are all above 30 and the per sequence GC content follows a normal distribution. Some of the sequences still have a high proportion of adapter content. 

### III) Trim 

Trim w/ [fastp](https://github.com/OpenGene/fastp)

```
nano fastp.sh

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="fastp_out_error"
#SBATCH --output="fastp_out"

echo "START"; date

module load fastp/0.19.7-foss-2018b 

base="/data/putnamlab/jillashey/BIO594_FinalProject"

for file in "$base"/data/raw/*.fastq
do
fastp -i $file -o $file.trim.fastp.fq 
done

mv *trim.fastp.fq "$base"/data/trimmed

echo "STOP"; date

sbatch fastp.sh
```

Submitted batch job 129676

### IV) Quality check trimmed reads 

Check number of files

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed
ls | wc -l
28
```

Count number of reads per file. Some reads have @HISEQ as header, some reads have @HWI

```
zgrep -c "HISEQ" *.trim.fastp.fq  > HISEQ_trim_length.txt

zgrep -c "HWI" *.trim.fastp.fq  > HWI_trim_length.txt
```

Run QC on raw reads 

```
nano fastqc_trim.sh 

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="fastqc_out_trim_error"
#SBATCH --output="fastqc_out_trim"

echo "START"; date

module load FastQC/0.11.9-Java-11 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

base="/data/putnamlab/jillashey/BIO594_FinalProject"

# Run fastqc on raw files 
for file in "$base"/data/trimmed/*.trim.fastp.fq 
do
fastqc $file --outdir "$base"/QC/trimmed/
done

# Compile with MultiQC
cd "$base"/QC/trimmed/
multiqc *fastqc.zip -o .

echo "STOP"; date

sbatch fastqc_trim.sh 
```

Submitted batch job 129499

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/trim_fastqc_per_base_sequence_quality_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/trim_fastqc_per_sequence_gc_content_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/trim_fastqc_adapter_content.png)

Like the raw QC, quality scores are all above 30 and the per sequence GC content still follows a normal distribution. But the adapters are gone!


### V) Align reads using different read aligners 

Obtain genomic/transcriptomic information for both species. 

##### A. cervicornis 
- Download genomic information [here](https://usegalaxy.org/u/skitch/h/acervicornis-genome)
- Path to genomic info: `/data/putnamlab/jillashey/genome/Acerv`

##### P. acuta 
- Download genomic information [here](http://cyanophora.rutgers.edu/Pocillopora_acuta/)
- Path to genomic info: `/data/putnamlab/jillashey/genome/Pacuta`

#### A) Align against genome - STAR

[STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference) is an aligner for RNA-seq data. It uses suffix arrays, seed clustering, and stitching. It can detect non-canonical splice sites, chimeric sequences, and also map full-length RNA sequences. It is fast, but memory intensive.

```
cd STAR 
```

First, use the genome to make an index

###### i) Acerv genome index

```
mkdir GenomeIndex_Acerv

nano GenomeIndex_Acerv.sh 

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="acerv_index_out_error"
#SBATCH --output="acerv_index_out"

echo "START"; date

module load STAR/2.7.2b-GCC-8.3.0  

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /data/putnamlab/jillashey/BIO594_FinalProject/STAR/GenomeIndex_Acerv \
--genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta \
--genomeSAindexNbases 13

echo "STOP"; date

sbatch GenomeIndex_Acerv.sh 
```

Submitted batch job 130259

###### ii) Pacuta genome index

```
mkdir GenomeIndex_Pacuta

nano GenomeIndex_Pacuta.sh 

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="pacuta_index_out_error"
#SBATCH --output="pacuta_index_out"

echo "START"; date

module load STAR/2.7.2b-GCC-8.3.0  

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /data/putnamlab/jillashey/BIO594_FinalProject/STAR/GenomeIndex_Pacuta \
--genomeFastaFiles /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.assembly.purged.fasta \
--genomeSAindexNbases 13

echo "STOP"; date

sbatch GenomeIndex_Pacuta.sh 
```

Submitted batch job 131712

`STAR` parameters: 

- `--runThreadN` - number of threads to run STAR
- `--runMode genomeGenerate` - runs a genome index generation job 
- `--genomeDir` - path to directory where genome indices are stored 
- `--genomeFastaFiles` - path to file with genome reference sequences 



Now align reads to genome 

###### iii) Acerv alignment 

```
mkdir AlignReads_acerv
cd AlignReads_acerv

nano AlignReads_Acerv.sh 

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="acerv_align_out_error"
#SBATCH --output="acerv_align_out"

echo "START"; date

module load STAR/2.7.2b-GCC-8.3.0 

F=/data/putnamlab/jillashey/BIO594_FinalProject/STAR/AlignReads_acerv 

# make symbolic link to trimmed data
ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .

# align reads 
array1=($(ls $F/*.trim.fastp.fq ))
for i in ${array1[@]}
do
STAR --runMode alignReads \
--genomeDir /data/putnamlab/jillashey/BIO594_FinalProject/STAR/GenomeIndex_Acerv/ \
--runThreadN 6 \
--readFilesIn ${i} \
--outFileNamePrefix ${i}. \
--outTmpDir ${i}_TMP \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outStd Log \
--twopassMode Basic \
--twopass1readsN -1 \
--outReadsUnmapped Fastx 
done

echo "STOP"; date

sbatch AlignReads_Acerv.sh 
```

Submitted batch job 131717

###### iv) Pacuta alignment 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/STAR
mkdir AlignReads_pacuta
cd AlignReads_pacuta

nano AlignReads_Pacuta.sh 

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="pacuta_align_out_error"
#SBATCH --output="pacuta_align_out"

echo "START"; date

module load STAR/2.7.2b-GCC-8.3.0 

F=/data/putnamlab/jillashey/BIO594_FinalProject/STAR/AlignReads_pacuta

# make symbolic link to trimmed data
ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .

# align reads 
array1=($(ls $F/*.trim.fastp.fq ))
for i in ${array1[@]}
do
STAR --runMode alignReads \
--genomeDir /data/putnamlab/jillashey/BIO594_FinalProject/STAR/GenomeIndex_Pacuta \
--runThreadN 6 \
--readFilesIn ${i} \
--outFileNamePrefix ${i}. \
--outTmpDir ${i}_TMP \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outStd Log \
--twopassMode Basic \
--twopass1readsN -1 \
--outReadsUnmapped Fastx 
done

echo "STOP"; date

sbatch AlignReads_Pacuta.sh 
```

Submitted batch job 131716 

`STAR` parameters: 

- `--runMode alignReads` - runs a read alignment job
- `--genomeDir` - path to genome index
- `--runThreadN` - number of threads to run STAR
- `--readFilesIn` - reading fastq files in 
- `--outFileNamePrefix` - naming the output files
- `--outTmpDir` - name of temporary directory that STAR creates (will be removed by STAR after its finished)
- `--outSAMtype` - type of SAM/BAM output (requesting unsorted and sorted BAM files). STAR is nice because it will just convert SAM to BAM files itself 
- `--outStd` - file type where output information will be stored
- `--twopassMode` - 2-pass mapping mode (Basic means all 1st pass junctions will be inserted into genome indices for 2nd pass)
- `--twopass1readsN` - number of reads to process for 1st pass (-1 means map all reads in 1st pass)
- `--outReadsUnmapped` - output of unmapped and partially mapped reads (Fastx indicates file type is a fasta/fastq file)


#### B) Align against transcriptome - Trinity & Bowtie2

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

First, create a de novo transcriptome w/ Trinity. Concatenate trimmed reads into single .fq file by species 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed

nano cat_all.sh

#!/bin/bash
#SBATCH --job-name="cat_all"
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="cat_all_out_error"
#SBATCH --output="cat_all_out"
#SBATCH -D /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed

echo "START"; date

cat 11_2.fastq.trim.fastp.fq 1_2.fastq.trim.fastp.fq 2_2.fastq.trim.fastp.fq 28_2.fastq.trim.fastp.fq 31_2.fastq.trim.fastp.fq 35_2.fastq.trim.fastp.fq 36_2.fastq.trim.fastp.fq 38_2.fastq.trim.fastp.fq 39_2.fastq.trim.fastp.fq 41_2.fastq.trim.fastp.fq 42_2.fastq.trim.fastp.fq 4_2.fastq.trim.fastp.fq 47_2.fastq.trim.fastp.fq > Pacuta_samples_all.fq

cat *Ac* > Acerv_samples_all.fq

echo "STOP"; date

sbatch cat_all.sh
```

Submitted batch job 130994

Move sample fq files to Trinity folder

```
mv Acerv_samples_all.fq Pacuta_samples_all.fq ../../Trinity/
```

Run Trinity for both species 

```
cd Trinity
```

###### i) Acerv Trinity 

```
nano trinity_acerv.sh

#!/bin/bash
#SBATCH --job-name="Trinity"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="trinity_acerv_out_error"
#SBATCH --output="trinity_acerv_out"
#SBATCH -D /data/putnamlab/jillashey/BIO594_FinalProject/Trinity
#SBATCH --mem=500GB

echo "START"; date

module load Trinity/2.9.1-foss-2019b-Python-3.7.4 

Trinity --seqType fq --single Acerv_samples_all.fq --SS_lib_type F --max_memory 125G --CPU 10 --output trinity_acerv --NO_SEQTK

echo "STOP"; date

sbatch trinity_acerv.sh
```

Submitted batch job 132511

###### ii) Pacuta Trinity 

```
nano trinity_pacuta.sh

#!/bin/bash
#SBATCH --job-name="Trinity"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="trinity_pacuta_out_error"
#SBATCH --output="trinity_pacuta_out"
#SBATCH -D /data/putnamlab/jillashey/BIO594_FinalProject/Trinity
#SBATCH --mem=500GB

echo "START"; date

module load Trinity/2.9.1-foss-2019b-Python-3.7.4 

Trinity --seqType fq --single Pacuta_samples_all.fq --SS_lib_type F --max_memory 125G --CPU 10 --output trinity_pacuta --NO_SEQTK

echo "STOP"; date

sbatch trinity_pacuta.sh
```

Submitted batch job 131000

`Trinity` parameters: 
- `--seqType` - fastq file type
- `--single` - single-end reads 
- `--SS_lib_type` - 
- `--max_memory` - 
- `--CPU` - 
- `--output` - name for output directory
- `--NO_SEQTK` - 

Nice, finished in about 4 days. Check out the `Trinity.fasta` file

```
cd trinity_pacuta

zgrep -c ">" Trinity.fasta
368329 # number of transcripts generated

head Trinity.fasta
>TRINITY_DN100155_c0_g1_i1 len=507 path=[0:0-506]
CTATACACAAACACCAATGATGAGGGATGGAGAATGATTCGTGTTGGAGATGACCCCCTAGGATTTGGTCAGTTTTGTCTTTTTTTTATTCTGATTTGAAATGAGCAAGCGAAAATTTGGTGAATAAGATGTTAGCTTAAGACTAGTGCATACAAACAGACACTCAACAGGCTGACACTGTTATGCATGTGGCTCAACCAGACATCAGACATTGTATAAGAAACGAGAAAATTGCTGATTTAGGTGACTGCCTAACTGATAACAGATAGGTGGATATAGTGGTGCACTGATTCACTGACTCAGACTGGCAGACTGACGGACTGATTGACTCAGACTGGCTGACTGACACACTGACCATCTGACCTACAGAATGACTGTCTGCATGACAAATTGATTGACTCACTGACTGATTCGCTCACCTATTCAGTCACTTGTTTTGCAAAGGTTATTAGAAAGGTGTGAACATTTGATACAAATCAGTCATAATAGAGGTGTGGTCTGCAAACC
>TRINITY_DN100130_c0_g1_i1 len=363 path=[0:0-362]
TTGTAGAACATGCTAGCCACCTCCATCATTTGCCCCATTACGTATCGCTGGAATCCTCGACGTTTAGATTTGTAGTGTAGCAGTAATCCCGTGCTTGCCTCTTCTGAACAGTAGAATGACGGTGACATGAGCTTTGGGTAAGCAAATCGCATGTGTTCGTGTAAATTGTCGATTCCTCGTAGGAAGTCACTGAAATGCCTGCCGCTTATTTGTATAAACCTATCGTAGCCGTAGTTGCTGAAAAATTTCACAAAGCACGTCCCAAAAAATTGAATGAAGTCTTCCGAGGTCATTTTTGTCTCATTTCCGAGTACTTCTGCCGCTGCAGACGCGATTTCCAACAACAGATTATCTGAGTACCGT
>TRINITY_DN100138_c0_g1_i1 len=259 path=[0:0-258]
CTCGGTGAAGTGCTCCGGATTCGGACCCTGGTCCCGGCAGTGCGGCTCTTGATGATTTCCGAAGAGGCGGCTTTGGGCCTGTCGCGAGTCGTGGGTGGTCTCCCATGTCATATTCGCACACGGCGTAAAGACATCCCACCCATGGGAGTGCAAACGGACGTTCTCGCAGGTGACTTGGCCGTGAAAATGCAAGGAATTCATGCAAGGGTCTGCTGGCACCTCTATCAAGATGTCAGAGGTGGAAAACGCACTGTGCGGA
>TRINITY_DN100156_c0_g1_i1 len=203 path=[0:0-202]
CAAATTGGTACCTCAATTTTGTTTCGAGCGTCGAAACGAGTTCTCCGGACCTCTCTCCGGATCATGTCTCCAATTTTCAGCCATGATCCCAATGTTGAGCCATGATCATTTTGAACATTTTCTTTCAAACATCTTTTTGCTCCCAATTGCCTGCGTAGCATGCTTTCAATTCGGAGAGCTATGTGCAAGAGCACATACCGTGG
>TRINITY_DN100195_c0_g1_i1 len=409 path=[0:0-408]
TTTTTTTTCCTCCACTTTTCTTTCCTTTTTTGCCTTTCTTTCCCTTTTTCTTTTTGTTGAGCTGTTTTCTGGCTCTATGCGCTCGCCATAGAGATTGTATTAATCTCGCCGCTCGTATTTTCAAATTCAGCTCTCTTTCAGCAGCTTCGGCCTTCTCTTTTGCAATTCTTCTCGCTTCCACGATTTCCAAGTACTCCTTCTCCAGGGTCTTGAAACGCTCCTCGAGTTCATTTAGCTGCCTTTTCTCCTCCGTGTAGATCGTGTCTAGGGCATCAAATTCATCTTGTCTTTCTCCCATATCTGTATCATACTTCTGTATCCAATTCTCAACTTCAGTTTCAATCTTGTACTTTCGCTTTCTCAAAGCCAGCTCACTCTCTCTATGAGTGAGTGTGTCGCTTTCAAGTTT
```
###### iii) Acerv Bowtie 

###### iv) Pacuta Bowtie 

Now run Bowtie w/ the `Trinity.fasta` as the Pacuta 'reference' transcriptome. First, build the bowtie index in the Bowtie2 directory

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Bowtie2
mkdir pacuta_index 

nano bowtie_index_pacuta.sh

#!/bin/bash
#SBATCH --job-name="Bowtie-index"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="bowtie_index_pacuta_out_error"
#SBATCH --output="bowtie_index_pacuta_out"
#SBATCH --mem=500GB

echo "START"; date

module load Bowtie2/2.4.4-GCC-11.2.0

bowtie2-build /data/putnamlab/jillashey/BIO594_FinalProject/Trinity/trinity_pacuta/Trinity.fasta pacuta_index/pacuta_idx

echo "STOP"; date

sbatch bowtie_index_pacuta.sh
```

Submitted batch job 132126

Once the index is made, use it to align reads to de novo transcriptome

```
mkdir pacuta_align

nano bowtie_align_pacuta.sh

#!/bin/bash
#SBATCH --job-name="Bowtie-align"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="bowtie_align_pacuta_out_error"
#SBATCH --output="bowtie_align_pacuta_out"
#SBATCH --mem=500GB

echo "START"; date

module load Bowtie2/2.4.4-GCC-11.2.0

F=/data/putnamlab/jillashey/BIO594_FinalProject/Bowtie2/pacuta_align

# make symbolic link to trimmed data
ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .

# align reads 
array1=($(ls $F/*.trim.fastp.fq ))
for i in ${array1[@]}
do
bowtie2 -q -x /data/putnamlab/jillashey/BIO594_FinalProject/Bowtie2/pacuta_index/pacuta_idx -U ${i} -S ${i}.sam
done

echo "STOP"; date

sbatch bowtie_align_pacuta.sh
```

Submitted batch job 132138

#### C) Pseudoalignment - Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/manual) is another tool to quantify RNA-seq data. Kallisto uses pseudoalignment to speed up the alignment process, meaning it can quantify reads without making actual alignments. It does this by identifying transcripts that a read is compatible with in order to quantify the transcript.

Build a kallisto index for both species 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

mkdir index
cd index

nano kallisto_index.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-index"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_index_out_error"
#SBATCH --output="kallisto_index_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

kallisto index -i acerv_index.idx /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.mRNA.fa 

kallisto index -i pacuta_index.idx /data/putnamlab/jillashey/genome/Pacuta/braker_v1.codingseq.fasta

echo "STOP"; date

sbatch kallisto_index.sh
```

Submitted batch job 131043

Use the index to pseudoalign the reads 

###### i) Acerv pseudoalignment 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto
mkdir Align_acerv
cd Align_acerv
```

Split samples up so that the samples that are 50 bp long are processed separately and the samples that are 120 bp long are also processed separately. This is important because Kallisto relies on fragment length to do its pseudoalignments.

| Sample Name         | Length | Species       |
| ------------------- | ------ | ------------- |
| 11\_2               | 120 bp | P.acuta       |
| 19\_T33\_Ac\_WK     | 50 bp  | A.cervicornis |
| 1\_2                | 50 bp  | P.acuta       |
| 24\_T12\_Ac\_FM     | 50 bp  | A.cervicornis |
| 25\_ctl1\_Ac\_GF\_2 | 123 bp | A.cervicornis |
| 27\_ctl2\_Ac\_YG\_2 | 123 bp | A.cervicornis |
| 28\_2               | 50 bp  | P.acuta       |
| 2\_2                | 50 bp  | P.acuta       |
| 31\_2               | 50 bp  | P.acuta       |
| 31\_T22\_Ac\_UV     | 50 bp  | A.cervicornis |
| 35\_2               | 50 bp  | P.acuta       |
| 35\_T43\_Ac\_MT     | 50 bp  | A.cervicornis |
| 36\_2               | 50 bp  | P.acuta       |
| 37\_T13\_Ac\_ML     | 50 bp  | A.cervicornis |
| 38\_2               | 50 bp  | P.acuta       |
| 38\_T23\_Ac\_IN     | 50 bp  | A.cervicornis |
| 39\_2               | 120 bp | P.acuta       |
| 41\_2               | 119 bp | P.acuta       |
| 41\_ctl3\_Ac\_RN\_2 | 123 bp | A.cervicornis |
| 42\_2               | 50 bp  | P.acuta       |
| 45\_T41\_Ac\_SC\_2  | 123 bp | A.cervicornis |
| 47\_2               | 50 bp  | P.acuta       |
| 47\_T31\_Ac\_JB     | 50 bp  | A.cervicornis |
| 4\_2                | 120 bp | P.acuta       |
| 52\_T11\_Ac\_II     | 50 bp  | A.cervicornis |
| 53\_T21\_Ac\_NH     | 50 bp  | A.cervicornis |
| 54\_T42\_Ac\_JQ     | 50 bp  | A.cervicornis |
| 57\_T32\_Ac\_NM     | 50 bp  | A.cervicornis |

Make symbolic links to trimmed data 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_acerv

ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .
```

Separate samples that have 50 or 120 bp into separate folders 

```
mkdir 50_bp 120_bp

# 50 bp
mv 28_2.fastq.trim.fastp.fq 38_2.fastq.trim.fastp.fq 1_2.fastq.trim.fastp.fq 31_2.fastq.trim.fastp.fq 38_T23_Ac_IN.fastq.trim.fastp.fq 47_2.fastq.trim.fastp.fq 19_T33_Ac_WK.fastq.trim.fastp.fq 31_T22_Ac_UV.fastq.trim.fastp.fq 47_T31_Ac_JB.fastq.trim.fastp.fq 35_2.fastq.trim.fastp.fq 2_2.fastq.trim.fastp.fq 52_T11_Ac_II.fastq.trim.fastp.fq 24_T12_Ac_FM.fastq.trim.fastp.fq 35_T43_Ac_MT.fastq.trim.fastp.fq 53_T21_Ac_NH.fastq.trim.fastp.fq 36_2.fastq.trim.fastp.fq 42_2.fastq.trim.fastp.fq 54_T42_Ac_JQ.fastq.trim.fastp.fq 37_T13_Ac_ML.fastq.trim.fastp.fq 57_T32_Ac_NM.fastq.trim.fastp.fq 50_bp

# 120 bp
mv 11_2.fastq.trim.fastp.fq 45_T41_Ac_SC_2.fastq.trim.fastp.fq 39_2.fastq.trim.fastp.fq 41_2.fastq.trim.fastp.fq 41_ctl3_Ac_RN_2.fastq.trim.fastp.fq 25_ctl1_Ac_GF_2.fastq.trim.fastp.fq 27_ctl2_Ac_YG_2.fastq.trim.fastp.fq 4_2.fastq.trim.fastp.fq 120_bp
```

Run kallisto quant step. 

For reads that are 50 bp: 

```
mkdir Align_acerv/output_50bp 
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_acerv/50_bp

nano kallisto_quant_50bp_acerv.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_50bp_acerv_out_error"
#SBATCH --output="kallisto_quant_50bp_acerv_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_acerv/50_bp/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/acerv_index.idx -o $F/Align_acerv/output_50bp --single -l 50 -s 5 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_50bp_acerv.sh 
```

Submitted batch job 131610

still not producing the .h5 file in output folder ?? Going to try with just one sample. still no...do i need this file anyway? just going to continue on. ah it seems that they may be phasing hdf5 out. see post [here](https://github.com/pachterlab/kallisto/releases)

For reads that are 120 bp: 

```
mkdir Align_acerv/output_120bp 
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_acerv/120_bp

nano kallisto_quant_120bp_acerv.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_120bp_acerv_out_error"
#SBATCH --output="kallisto_quant_120bp_acerv_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_acerv/120_bp/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/acerv_index.idx -o $F/Align_acerv/output_120bp --single -l 120 -s 10 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_120bp_acerv.sh 
```

Submitted batch job 131611

###### ii) Pacuta pseudoalignment 

Make directories and create symbolic links to samples 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto
mkdir Align_pacuta
cd Align_pacuta
mkdir 50_bp 120_bp

ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .
```

Separate samples that have 50 or 120 bp into separate folders 

```
cd 50_bp

# 50 bp
mv 28_2.fastq.trim.fastp.fq 38_2.fastq.trim.fastp.fq 1_2.fastq.trim.fastp.fq 31_2.fastq.trim.fastp.fq 38_T23_Ac_IN.fastq.trim.fastp.fq 47_2.fastq.trim.fastp.fq 19_T33_Ac_WK.fastq.trim.fastp.fq 31_T22_Ac_UV.fastq.trim.fastp.fq 47_T31_Ac_JB.fastq.trim.fastp.fq 35_2.fastq.trim.fastp.fq 2_2.fastq.trim.fastp.fq 52_T11_Ac_II.fastq.trim.fastp.fq 24_T12_Ac_FM.fastq.trim.fastp.fq 35_T43_Ac_MT.fastq.trim.fastp.fq 53_T21_Ac_NH.fastq.trim.fastp.fq 36_2.fastq.trim.fastp.fq 42_2.fastq.trim.fastp.fq 54_T42_Ac_JQ.fastq.trim.fastp.fq 37_T13_Ac_ML.fastq.trim.fastp.fq 57_T32_Ac_NM.fastq.trim.fastp.fq 50_bp

# 120 bp
mv 11_2.fastq.trim.fastp.fq 45_T41_Ac_SC_2.fastq.trim.fastp.fq 39_2.fastq.trim.fastp.fq 41_2.fastq.trim.fastp.fq 41_ctl3_Ac_RN_2.fastq.trim.fastp.fq 25_ctl1_Ac_GF_2.fastq.trim.fastp.fq 27_ctl2_Ac_YG_2.fastq.trim.fastp.fq 4_2.fastq.trim.fastp.fq 120_bp
```

now can run kallisto quant step. 

For reads that are 50 bp: 

```
mkdir Align_pacuta/output_50bp 
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_pacuta/50_bp

nano kallisto_quant_50bp_pacuta.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_50bp_pacuta_out_error"
#SBATCH --output="kallisto_quant_50bp_pacuta_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_pacuta/50_bp/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/pacuta_index.idx -o $F/Align_pacuta/output_50bp --single -l 50 -s 5 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_50bp_pacuta.sh 
```

Submitted batch job 131612

```
mkdir Align_pacuta/output_120bp 
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_pacuta/120_bp

nano kallisto_quant_120bp_pacuta.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_120bp_pacuta_out_error"
#SBATCH --output="kallisto_quant_120bp_pacuta_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_pacuta/120_bp/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/pacuta_index.idx -o $F/Align_pacuta/output_120bp --single -l 120 -s 10 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_120bp_pacuta.sh 
```

Submitted batch job 131613

Pseudoalignment information by species:

| Species Alignment | k-mer length | # of targets       | # of k-mers       | # of equivalence classes |
| --------------| ------ | -------- | ------ | -------- | -------- |
| A.cervicornis | 31 | 33,322 | 34,006,771 | 84,017 |
| P.acuta | 31 | 38,913 | 38,267,300 | 88,498 |

### VI) Compare 

Keep in mind, none of these can be DIRECTLY compared, as they are all using different algorithms to calculate alignments. But it is still interesting to qualatitvely compare across aligner types. 

Sample 31_2 (P.acuta) had a very low # of read counts (~6000 raw reads), so it will be excluded from comparisons. Additionally, sample 24_T12_Ac_FM (A.cervicornis) will be excluded because it had strangely low alignment for an Acerv sample. I'm thinking that it was mislabeled, as these samples were processed with several other species. Both 31_2 and 24_T12_Ac_FM are highlighted in yellow on the alignment csv files. 

#### Average alignment rates for each tool 

| Species | STAR | Bowtie | Kallisto |
| --------------| ------ | -------- | ------ | -------- | 
| A.cervicornis | 70.72 | xx | 34.67 | 
| P.acuta | 69.40 | 94.93 | 33.54 | 

These are interesting results. I was not expecting the Bowtie alignment to be so high, especially when I used a de novo transcriptome. Because I used a de novo transcriptome, that does give me a little less confidence in the alignment quality. 

See STAR alignment csv, Bowtie alignment csv and Kallisto alignment csv for more detail. 

#### Time, memory used, and CPU efficiency for each tool

I am also interested in how much time it take for each tool to process my samples, as well as how much memory and CPU each tool requires. For these calculations, I am factoring in how long it took to build the indices, as well as align. I had 28 total samples that were between 1.5G to 5.5G in size. 

##### STAR

Total Acerv time: 2:33:05
Total Acerv memory utilized: 15.91 GB
Total Acerv CPU utilized: 13:16:04

Total Pacuta time: 4:04:38
Total Pacuta memory utilized: 15.6 GB
Total Pacuta CPU utilized: 22:27:01

##### Bowtie + Trinity 

Total Acerv time: xx
Total Acerv memory utilized: xx
Total Acerv CPU utilized: xx

Total Pacuta time: 4-13:30:17
Total Pacuta memory utilized: 124:45 GB
Total Pacuta CPU utilized: 10-08:51:54

##### Kallisto

Total Acerv time: 0:47:42
Total Acerv memory utilized: 4.85 GB
Total Acerv CPU utilized: 0:44:22

Total Pacuta time: 00:47:18
Total Pacuta memory utilized: 4.89 GB
Total Pacuta CPU utilized: 0:45:18


### VII) References 

- [Bray et al. 2016](https://www.nature.com/articles/nbt.3519) - Near-optimal probabilistic RNA-seq quantification
- [Dobin et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) - STAR: ultrafast universal RNA-seq aligner
- [Haas et al. 2013](https://www.nature.com/articles/nprot.2013.084) - De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis
- [Langmead et al. 2009](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25) - Ultrafast and memory-efficient alignment of short DNA sequences to the human genome