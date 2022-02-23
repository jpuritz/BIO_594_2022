# BIO 594 Exercise Week 4
## `fastp` creating and filtering genomics data sets

GitHub documentation for `fastp` can be found [here](https://github.com/OpenGene/fastp#features)

[Documentation](http://bio-bwa.sourceforge.net/bwa.shtml) for bwa alignment package used in the Map.Ex activity in BIO 594

### Goal from BIO 594 Github for this markdown
> Create filtered read sets for all data in /home/BIO594/DATA/Week4/realdata/
>  * Post the filtered read sets in a directory within the class github repository /Exercises/Week04/ directory

#### Start -Data into working directory
Once in Kitt enter the wk4 directory and make a new directory to practice in, then make a data directory within it.
````
```bash
cd wk4
mkdir goal
mkdir data
```
````
Inspect directory with data then copy data to the "goal" directory
````
```bash
cd /home/BIO594/DATA/Week4/realdata/
cp -R /home/BIO594/DATA/Week4/realdata/* ~/wk4/goal/data # -R flag for copying a directory, * wildcard specifying all files
```
````
There are four datasets: Exome_capture, rad_seq, rna_seq, and wgs

#### Loading in fastp

Check to see if fastp is already installed to conda in kitt

````
```bash
conda list
```
````
output - fastp is loaded already on kitt
```
...
fastp                     0.12.4                        0    bioconda
...
```
Create a conda envionment with fastp and activate it
````
```bash
conda info --envs #checks currently available environments
conda create -n fast-env fastp  # flag "environment name" "names of packages to load
conda activate fastp-env # starts new environment
```
````


#### Attmepting fastp directly on data file – no conversion to SAM/BAM format
Using default quality filter settings for fastp – no additional flags beyond input and output

##### Exome capture dataset
````
```bash
mkdir ~/wk4/goal/output/exom_capture #making output directory
cd data/exome_capture
fastp -i Capture1.F.fq -o ../../output/exome_capture/cp1.f.fq # -i "input file" -o "output file name"
```
````
Files are output as named to specified folder. The fastp.html and fastp.json files...reports(?) are left in the folder where the data came was stored `~/data/wk4/goal/data/exome_capture`. Opening the output file as is or adding the .html extension does not produce the report on the [fastp documentation page](https://github.com/OpenGene/fastp#features). Otherwise it seems to have worked. 

***Update*** the reports `fastp.html` and `fastp.json` open when they are selected from a program like `WinSCP` and specified to open with a web browser. 
* Ask how to specify where these files output and how to name them

In the mean time I will manually move and rename these manually using `WinSCP`

````
```bash
cd data/exome_capture
fastp -i Capture1.R.fq -o ../../output/exome_capture/cp1.R.fq # -i "input file" -o "output file name"
```
````
Paired end `fastp` with gzipped files
````
```bash
fastp -i CASE_J03.F.fq.gz -I CASE_J03.R.fq.gz  -o ../../output/exome_capture/Out_CASE_J03.F.fq.gz  -O ../../output/exome_capture/Out_CASE_J03.R.fq.gz
```
````

#### rad_seq datasets

````
```bash
cd ../rad_seq

fastp -i FCLib1.F.fastq.gz -I FCLib1.R.fastq.gz  -o ../../output/rad_seq/Out_FCLib1.F.fastq.gz  -O ../../output/rad_seq/Out_FCLib1.R.fastq.gz

ls
```
````
Move and rename report files

````
```bash

fastp -i WC1_407.F.fq.gz -I WC1_407.R.fq.gz  -o ../../output/rad_seq/Out_WC1_407.F.fq.gz  -O ../../output/rad_seq/Out_WC1_407.R.fq.gz

ls
```
````
#### rna_seq data
Adding the -g flag to specify trimming of poly G tails. `fastp` documentation says this happens with rna_seq data. *This may be redundant* since poly G trimming is eneabled by default.

````
```bash

cd ../rna_seq

fastp -g -i rna1.F.fq.gz -I rna1.R.fq.gz  -o ../../output/rna_seq/Out_rna1.F.fq.gz  -O ../../output/rna_seq/Out_rna1.R.fq.gz

ls
```
````
Move and rename report files to output file

#### wgs datasets
````
```bash
cd ../wgs

fastp -i HC_1.F.fq.gz -I HC_1.R.fq.gz  -o ../../output/wgs/Out_HC_1.F.fq.gz  -O ../../output/wgs/Out_HC_1.R.fq.gz

ls
```
````
Move and rename report files

Work complete put in specified folder, pull then push to git

Read to understand all the parameters you are messing with.