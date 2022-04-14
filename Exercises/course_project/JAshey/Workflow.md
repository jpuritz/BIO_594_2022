Code for my bio594 project 

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
mkdir data QC STAR Bowtie2 Kallisto
cd data
mkdir raw trimmed
```

Make links to Acerv and Pacuta data 

```
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

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/fastqc_per_base_sequence_quality_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/fastqc_per_sequence_gc_content_plot.png)

![](https://raw.githubusercontent.com/jpuritz/BIO_594_2022/main/Exercises/course_project/JAshey/images/fastqc_adapter_content_plot.png)

I didn't add all the plots, just the ones that I thought were most important. Sequences were either 50 or 125 bp long. Quality scores are all above 30 and the per sequence GC content follows a normal distribution. Some of the sequences still have a high proportion of adapter sequences. 

### Trim 

Trimmomatic will be used to trim and clean reads. 

Copy Ilumina adapter clip information into directory. This info is from Ilumina about the adapter lengths that were specifically used to sequence these samples. 

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/data/raw
cp Illumina_adapter_reads_PE_SE.fa /data/putnamlab/jillashey/BIO594_FinalProject/data
cd /data/putnamlab/jillashey/BIO594_FinalProject/data/raw
```


redo trimming step, check adapter code 4/13/22

try trimming w/ [fastp](https://github.com/OpenGene/fastp)

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

Rerun QC with trimmed (fastp) data














Run Trimmomatic 

```
nano trimmomatic.sh

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="trimmomatic_out_error"
#SBATCH --output="trimmomatic_out_error"

echo "START"; date

module load Trimmomatic/0.39-Java-11 

base="/data/putnamlab/jillashey/BIO594_FinalProject"

for file in "$base"/data/raw/*.fastq
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -phred33 $file $file.trim.fq ILLUMINACLIP:"$base"/data/Illumina_adapter_reads_PE_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> amtTrimmed.txt
done 

mv *trim.fq "$base"/data/trimmed

echo "STOP"; date

sbatch trimmomatic.sh
```

Submitted batch job 129350

`Trimmomatic` arguments:

- `SE` = single end reads 
- `phred33` = specifies base quality encoding 
- `ILLUMINACLIP` = cuts adapters and other Illumina specific sequences from read 
- `LEADING` = cut off bases (in this case, up to 3) at start of read, if below a certain quality 
- `TRAILING` = cut off bases (in this case, up to 3) at end of read, if below a certain quality 
- `SLIDINGWINDOW` = xxxxxxx
- `MINLEN` = drop read if below a certain length (in this case, 20 bp)

### Quality check trimmed reads 

Check number of files to make sure everything was moved properly

```
cd /data/putnamlab/jillashey/BIO594_FinalProject/data/trimmed
ls | wc -l
28
```

Count number of reads per file. Some reads have @HISEQ as header, some reads have @HWI

```
zgrep -c "HISEQ" *.trim.fq > HISEQ_trim_length.txt

zgrep -c "HWI" *.trim.fq > HWI_trim_length.txt
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
for file in "$base"/data/trimmed/*.trim.fq
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

Copy MultiQC files onto computer and look at plots

### Align reads using different read aligners 

Obtain genomic/transcriptomic information for both species. 

##### A. cervicornis 
- Download genomic information here xxxx
- Path to genomic info: `/data/putnamlab/jillashey/genome/Acerv`

##### P. acuta 
- Download genomic information here xxxx
- Path to genomic info: `/data/putnamlab/jillashey/genome/Pacuta`

#### Align against genome - STAR

```
cd STAR 
```

##### A. cervicornis 

Use the genome to make an index

```
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

##### P. acuta 

Use the genome to make an index

```
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
--genomeFastaFiles /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.assembly.purged.fasta

echo "STOP"; date

sbatch GenomeIndex_Pacuta.sh 
```


#### Align against transcriptome - Bowtie2

#### Pseudoalignment - Kallisto

### Compare 