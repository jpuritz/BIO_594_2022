Code for my bio594 project. Analysis done on Andromeda 

Proposal [here](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/course_project/JAshey/Ashey_BIO594_FinalProjectProposal.md)

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
--genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta

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

Submitted batch job 130432

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
--outFileNamePrefix ${i}.
--outTmpDir ${i}_TMP
--outSAMtype BAM Unsorted SortedByCoordinate \
--outStd Log BAM_Unsorted BAM_Quant \
--quantMode TranscriptomeSAM \
--twopassMode Basic \
--twopass1readsN -1 \
--outReadsUnmapped Fastx 
done

echo "STOP"; date

sbatch AlignReads_Acerv.sh 
```

Submitted batch job 130433

Convert SAM to BAM files using samtools 

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
--genomeDir /data/putnamlab/jillashey/BIO594_FinalProject/STAR/GenomeIndex_Pacuta/ \
--runThreadN 6 \
--readFilesIn ${i} \
--outFileNamePrefix ${i}.
--outTmpDir ${i}_TMP
--outSAMtype BAM Unsorted SortedByCoordinate \
--outStd Log BAM_Unsorted BAM_Quant \
--quantMode TranscriptomeSAM \
--twopassMode Basic \
--twopass1readsN -1 \
--outReadsUnmapped Fastx 
done

echo "STOP"; date

sbatch AlignReads_Pacuta.sh 
```

Submitted batch job 131040

`STAR` parameters: 

- `runMode alignReads` - runs a read alignment job
- `genomeDir` - path to genome index
- `runThreadN` - number of threads to run STAR
- `readFilesIn` - reading fastq files in 
- `outFileNamePrefix` - naming the output files
- `outTmpDir` - name of temporary directory STAR creates to align reads. This will be deleted by the end of the job
- `outSAMtype`
- `outStd`
- `quantMode`
- `twopassMode`
- `twopass1readsN`
- `outReadsUnmapped`

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
- `seqType` - 
- `single` - 
- `SS_lib_type` - 
- `max_memory` - 
- `CPU` - 
- `output` - 



```
## Kevin script
Trinity --cite
Trinity --version
Trinity --seqType fq --max_memory 125G --bflyCalculateCPU --CPU 10 --FORCE --left ../raw_reads/Sample_R1_all.fastq.gz --right ../raw_reads/Sample_R2_all.fastq.gz --SS_lib_type RF --full_cleanup #run assembly --output ./Pastreoides_holobiont_assembly
```

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

###### Pacuta pseudoalignment 


### Compare 