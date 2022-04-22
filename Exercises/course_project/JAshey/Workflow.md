Code for my bio594 project. Analysis done on Andromeda 

Proposal [here](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Ashey_BIO594_FinalProjectProposal.md)

idk who said this but wise words: "the computer is never wrong, it is you who is wrong"

### Set up directories and data 

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
zgrep -c "HISEQ" *.fastq > HISEQ_raw_length.txt

zgrep -c "HWI" *.fastq > HWI_raw_length.txt
```

### Quality check

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

### Trim 

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

### Quality check trimmed reads 

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


### Align reads using different read aligners 

Obtain genomic/transcriptomic information for both species. 

##### A. cervicornis 
- Download genomic information [here](https://usegalaxy.org/u/skitch/h/acervicornis-genome)
- Path to genomic info: `/data/putnamlab/jillashey/genome/Acerv`

##### P. acuta 
- Download genomic information [here](http://cyanophora.rutgers.edu/Pocillopora_acuta/)
- Path to genomic info: `/data/putnamlab/jillashey/genome/Pacuta`

#### Align against genome - STAR

[STAR](https://github.com/alexdobin/STAR)

```
cd STAR 
```

##### Use the genome to make an index

###### Acerv genome index

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

###### Pacuta genome index

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

- `runThreadN` - number of threads to run STAR
- `runMode genomeGenerate` - runs a genome index generation job 
- `genomeDir` - path to directory where genome indices are stored 
- `genomeFastaFiles` - path to file with genome reference sequences 



##### Align reads to genome 

###### Acerv alignment 

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

Submitted batch job 131682

rerunning star 4/21/22


Didn't make bam files???? so going to convert SAM to BAM files using [samtools](https://www.htslib.org/doc/samtools-sort.html)
```
nano samtools.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="samtools_out_error"
#SBATCH --output="samtools_out"

echo "START"; date

module load SAMtools/1.9-foss-2018b

F=/data/putnamlab/jillashey/BIO594_FinalProject/STAR/AlignReads_acerv

array1=($(ls $F/*Aligned.out.sam))
for i in ${array1[@]}
do
samtools sort -@ 8 -o ${i}.bam ${i}
done

echo "STOP"; date

sbatch samtools.sh
```

Submitted batch job 131037

###### Pacuta alignment 

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

Submitted batch job 131717 -- rerunning star

tried: rerunning genome index, now taking out --outSAMtype


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

convert SAM to BAM files using [samtools](https://www.htslib.org/doc/samtools-sort.html)

```
nano samtools.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="samtools_out_error"
#SBATCH --output="samtools_out"

echo "START"; date

module load SAMtools/1.9-foss-2018b

F=/data/putnamlab/jillashey/BIO594_FinalProject/STAR/AlignReads_pacuta

array1=($(ls $F/*Aligned.out.sam))
for i in ${array1[@]}
do
samtools sort -@ 8 -o ${i}.bam ${i}
done

echo "STOP"; date

sbatch samtools.sh
```

Submitted batch job 131285

`samtools` parameters: 

- `sort` - sort alignment by coordinate and create bam file
- `-@` - number of sorting and compression threads 
- `-o` - name of output file


#### Align against transcriptome - Trinity & Bowtie2

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

###### Acerv Trinity 

###### Pacuta Trinity 

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

Trinity --seqType fq --single Pacuta_samples_all.fq --SS_lib_type F --max_memory 125G --CPU 10 --output trinity_pacuta

echo "STOP"; date

sbatch trinity_pacuta.sh
```

Submitted batch job 131000

`Trinity` parameters: 
- `seqType` - fastq file type
- `single` - single-end reads 
- `SS_lib_type` - 
- `max_memory` - 
- `CPU` - 
- `output` - name for output directory

#### Pseudoalignment - Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/manual)

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

###### Acerv pseudoalignment 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto
mkdir Align_acerv
cd Align_acerv
mkdir quant_acerv

nano kallisto_quant_acerv.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_acerv_out_error"
#SBATCH --output="kallisto_quant_acerv_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# make symbolic link to trimmed acerv data
ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*Ac*trim.fastp.fq .

# align reads 
array1=($(ls $F/Align_acerv/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/acerv_index.idx -o $F/Align_acerv/quant_acerv --single -l 85 -s 40 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_acerv.sh 
```

Submitted batch job 131288. 

Some reads are 50 bp, some are ~120. That is why the length is 85 and the std deviation is 40

hmm i got a Kallisto output, but the .hd5 file is not there. The alignment rate also appears to be pretty low and there are a lot of contigs, k-mers, and equivalence classes. 

Split samples up so that the samples that are 50 bp long are processed separately and the samples that are 120 bp long are also processed separately. 

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

now can run kallisto quant step. 

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

still not producing the .h5 file in output folder ?? Going to try with just one sample 

```
kallisto quant -i /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/index/acerv_index.idx -o /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_acerv/test --single -l 50 -s 5 --plaintext 19_T33_Ac_WK.fastq.trim.fastp.fq
```

still no...do i need this file anyway? just going to continue on. ah it seems that they may be phasing hdf5 out. see post [here](https://github.com/pachterlab/kallisto/releases)

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

###### Pacuta pseudoalignment 

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

Pseudoalignment by species:

| Species Alignment | k-mer length | # of targets       | # of k-mers       | # of equivalence classes |
| --------------| ------ | -------- | ------ | -------- | -------- |
| A.cervicornis | 31 | 33,322 | 34,006,771 | 84,017 |
| P.acuta | 31 | 38,913 | 38,267,300 | 88,498 |

The alignment values are pretty low in both species. maybe just run with only one species? ie Acerv run with only acerv samples and index, Pacuta run with only pacuta samples and index

Make directories, make symbolic link to samples, put samples in proper directories 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto
mkdir Align_acerv/acerv_only Align_pacuta/pacuta_only

ln -s /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed/*trim.fastp.fq .

# Acerv - move first 
mv *Ac* Align_acerv/acerv_only

# Pacuta
mv *.trim.fastp.fq Align_pacuta/pacuta_only
```

Run Kallisto for only Acerv samples 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_acerv
mkdir output_acerv_only
cd acerv_only

nano kallisto_quant_acerv_only.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_acerv_only_out_error"
#SBATCH --output="kallisto_quant_acerv_only_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_acerv/acerv_only/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/acerv_index.idx -o $F/Align_acerv/output_acerv_only --single -l 85 -s 40 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_acerv_only.sh 
```

Going back to the weird 85 fragment length bc the fragments are either 120 or 50 bp. Submitted batch job 131620

Run Kallisto for only Pacuta samples 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/Kallisto/Align_pacuta
mkdir output_pacuta_only
cd pacuta_only

nano kallisto_quant_pacuta_only.sh

#!/bin/bash
#SBATCH --job-name="Kallisto-quant"
#SBATCH -t 336:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="kallisto_quant_pacuta_only_out_error"
#SBATCH --output="kallisto_quant_pacuta_only_out"
#SBATCH --mem=500GB

echo "START"; date

module load kallisto/0.46.2-foss-2020b

F=/data/putnamlab/jillashey/BIO594_FinalProject/Kallisto

# align reads 
array1=($(ls $F/Align_pacuta/pacuta_only/*.trim.fastp.fq ))
for i in ${array1[@]}
do
kallisto quant -i $F/index/pacuta_index.idx -o $F/Align_pacuta/output_pacuta_only --single -l 85 -s 40 -b 30 ${i}
done

echo "STOP"; date

sbatch kallisto_quant_pacuta_only.sh 
```

Going back to the weird 85 fragment length bc the fragments are either 120 or 50 bp. Submitted batch job 131621

hmm interesting. When running species separately, there are less reads that are pseudoaligned in both species. I'm going to stick with my 50bp and 120bp analysis because Kallisto is relying on fragment length

### Compare 

Keep in mind, none of these can be DIRECTLY compared

comparing kallisto and star github [post](https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md)