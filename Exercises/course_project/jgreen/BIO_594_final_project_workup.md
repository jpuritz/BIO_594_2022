# BIO 594 final project workup

## Problem Statement

EecSeq is a new tool for evolutionary biologists to explore questions concerning genes under selection. By reducing the genomic information captured we can greatly expand the number of samples a sequencing study design can accommodate. Additionally EecSeq has no need for a reference to design probes opening up the potential to study non-model or obscure model organisms, but without a reference EecSeq reads must be assembled through a de novo bioinformatic framework. The EecSeq protocol creates cDNA libraries of high quality that could be sequenced and guide EecSeq assembly through a de novo transcriptome using these previously unsequenced cDNA libraries. This project will explore de novo assembly of RNAseq reads into an annotated de novo transcriptome that will act as a sudo-reference for EecSeq read assembly.

## Goals

1) Link and organize reads
2) Process reads
       * Trim adapters
       * Assess trimming
       * Unique read filtering or....
       * Nommalization
3) Cluster cDNA files
4) Use itero to assemble gDNA reads with cDNA "seeds"
5) Use trinity to build a de novo transcriptome
6) Assess the assembly
       * BUSCO
       * Quast
       * Detonate
       * Transrate
       * Trinity Stats

Using https://itero.readthedocs.io/en/latest/ to assembly gDNA reads. A simple PE read stitch for the cDNA (something like https://cme.h-its.org/exelixis/web/software/pear/). Use this uniquely filtered cDAN file as a single set of unique merges as "seeds" for itero.  

### Sources of Data

* cDNA reads from CASE study
* gDNA reads from CASE study


## cDNA reads

Make directory for cDNA reads

```bash
mkdir cDNA
```

Symbolic link for reads in the cDNA folder for cDNA reads

```bash
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool1/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool2/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool3/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool4/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool4/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool5/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool6/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool7/
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/cDNA/pool8/
```

## Process reads

### Fastp filtering

fastp filtering

```bash
# Options for fastp used are
# --thread: number of threads
# -V: verbose log information
# -e: average quality and if anything is below that value remove it
# -q: quality value of that defined base
# -p: overrepresentation analysis
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p -i /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1/cDNApool-1_R1_001.fastq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1/cDNApool-1_R2_001.fastq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1_trimmed/cDNApool-1_R1_001_trim.fastq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1_trimmed/cDNApool-1_R2_001.trim.fastq.gz
```

Run fastqc on all cDNA trimmed reads

```bash
fastqc -f fastq *.fq.gz
multiqc -o multiqc .
```

Upon doing QC for the raw samples we found an adundance of adapters that fastp was not recognizing even if we would feed fastp a .fa file with those specific adapters added to it. So I will try demultiplexing these samples

Need to demultiplex cDNA samples. I'll work on demultiplexing just one pool for right now

```bash
process_shortreads -1 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1/cDNApool-1_R1_001.fastq.gz -2 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool1/cDNApool-1_R2_001.fastq.gz -b /home/jgreen/repos/BIO594_work/course_project/data/cDNA/barcodes.txt --inline_inline -D --barcode_dist_1 2 --barcode_dist_2 2 -o /home/jgreen/repos/BIO594_work/course_project/data/cDNA/demux/pool1/
```

We can create a for loop to demultiplex all of our samples 

```bash
for ((i=1;i<=8;i++));
do
process_shortreads -1 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool$i/cDNApool-${i}_R1_001.fastq.gz -2 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool$i/cDNApool-${i}_R2_001.fastq.gz -b /home/jgreen/repos/BIO594_work/course_project/data/cDNA/barcodes.txt --inline_inline -D --barcode_dist_1 2 --barcode_dist_2 2 -o /home/jgreen/repos/BIO594_work/course_project/data/cDNA/demux/pool$i/
done
```

Alright we are abandoning that and just linking the demultiplexed samples. Should've done that from the beginning

Symbolic link for demultiplexed cDNA reads

```bash
ln -s /RAID_STORAGE2/Shared_Data/CASE_RNA/demux/
```

### Fastp demultiplexed reads

Use fastp to trim the sequences. We are using a neat trick I found where we can declare string arrays and use for loops to go through that list. This is especially helpful when we have a ton of smaples.

```bash
# Set string array with the prefix name of each read file.
declare -a StringArray4=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")

# For loop to print the string array just to confirm the proper names
for i in "${StringArray4[@]}"
do
echo $i
done

# For loop that loops through the string array names and sets is as the variable i. To call the variable use the ${i}. 
# Options for fastp used are
# --thread: number of threads
# -V: verbose log information
# -e: average quality and if anything is below that value remove it
# -q: quality value of that defined base
# -p: overrepresentation analysis
for i in "${StringArray[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/cDNA/demux/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/cDNA/demux/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/${i}_trim.R.fq.gz
done
```

```bash
fastqc -f fastq *.fq.gz
multiqc -o multiqc .
```

We are going to explore some ways to reduce the data going into the clustering. The first method will be merging the trimmed reads using PEAR and clustering the merged reads to reduce redundancy across experiments. The second method will be uniq read filtering which will give us one file to cluster and pipe into cdhit. The third method will be using in silico read normalization script in trinity to process the trimmed reads and then use PEAR to merge the reads. After merging we will use a experiment specific design to create clustered reads to a set as seeds for itero.

### Method 1: Merging and Clustering

From here we will be merging the forward and reverse reads

```bash
# use conda to create pear env
# mamba create -n pear pear
# conda activate pear
base="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/merge/"
declare -a StringArray=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
for i in "${StringArray[@]}"
do
pear -j 20 -f "$i"_trim.F.fq.gz -r "$i"_trim.R.fq.gz -o "$base""$i"_trim_merge.fq.gz
done
```

Clustering reads

```bash
# use conda to create cdhit env
# mamba create -n cdhit cd-hit
# conda activate cdhit
declare -a StringArray=("CASE" "CA" "CON" "SE")
cd-hit-est -i "$i".cDNA.assembled.fastq -o "$i".cDNA.c95.cluster -c 0.95 -n 10 -d 0 -M 16000 -T 8
```

### Method 2: Unique read filtering

This series of commands are taken from the ddocent tutorial on unique read filtering

```bash
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"

cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile

more uniqseq.data

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [2:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF
```

                       Number of Unique Sequences with More than X Coverage (Counted within individuals)

    4e+07 ++----------+-----------+----------+-----------+-----------+-----------+----------+-----------+----------++
          *           +           +          +           +           +           +          +           +           +
          *                                                                                                         |
  3.5e+07 +*                                                                                                       ++
          |*                                                                                                        |
          | *                                                                                                       |
    3e+07 ++*                                                                                                      ++
          | *                                                                                                       |
          |  *                                                                                                      |
  2.5e+07 ++ *                                                                                                     ++
          |  *                                                                                                      |
    2e+07 ++  *                                                                                                    ++
          |   *                                                                                                     |
          |    *                                                                                                    |
  1.5e+07 ++   *                                                                                                   ++
          |    *                                                                                                    |
          |     *                                                                                                   |
    1e+07 ++    ***                                                                                                ++
          |        **                                                                                               |
          |          *                                                                                              |
    5e+06 ++          ******                                                                                       ++
          |                 ******                                                                                  |
          +           +           *****************************************************     +           +           +
        0 ++----------+-----------+----------+-----------+-----------+-----------+-----******************************
          2           4           6          8           10          12          14         16          18          20
                                                           Coverage


```bash
parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv

for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
rm ufile

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF
```

  900000 ++-----------+-------------+------------+-------------+------------+------------+-------------+-----------++
         +            +             +            +             +            +            +             +            +
         ****                                                                                                       |
  800000 ++  **                                                                                                    ++
         |     **                                                                                                   |
  700000 ++      **                                                                                                ++
         |         **                                                                                               |
         |           *                                                                                              |
  600000 ++           ******                                                                                       ++
         |                  ***                                                                                     |
  500000 ++                    ****                                                                                ++
         |                         *                                                                                |
         |                          **********                                                                      |
  400000 ++                                   ***                                                                  ++
         |                                       *******                                                            |
  300000 ++                                             *****                                                      ++
         |                                                   **                                                     |
         |                                                     **********                                           |
  200000 ++                                                              ***                                       ++
         |                                                                  *************                           |
  100000 ++                                                                              **************            ++
         |                                                                                             **************
         +            +             +            +             +            +            +             +            +
       0 ++-----------+-------------+------------+-------------+------------+------------+-------------+-----------++
         2            3             4            5             6            7            8             9            10
                                                    Number of Individuals
```

```bash
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs
```

439323 uniq.k.4.c.4.seqs

Convert into Fasta file

```bash
cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
```

Now we will extract the forward reads and cluster them

```bash
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1
``` 

When we take a look at the number of sequences retained from this process we have a total of **4,234 cDNA sequence**

This is not nearly enough to fully assemble the exome. We will be moving on from this method to in silico read normalization.

### Method 3: Insilico read normalization, merging, clustering

```bash
# use conda to create trinity env
# mamba create -n trinity trinity
# conda activate trinity
base="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/norm/"
reads="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/"
declare -a StringArray=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
for i in "${StringArray[@]}"
do
insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 30 --left "$i"_trim.F.fq.gz --right "$i"_trim.R.fq.gz --pairs_together --PARALLEL_STATS --CPU 10 --output "$base"
done
```

```bash
for fq in *.fq
do
echo $fq
cat $fq | echo $((`wc -l`/4))
done
```
|Sample Name|Read Count|
|-----------|----------|
|CA_J06_trim_norm.F|5003851|
|CA_J06_trim_norm.R|5003851|
|CA_J08_trim_norm.F|5514394|
|CA_J08_trim_norm.R|5514394|
|CA_J11_trim_norm.F|6112644|
|CA_J11_trim_norm.R|6112644|
|CA_J18_trim_norm.F|7682106|
|CA_J18_trim_norm.R|7682106|
|CASE_J03_trim_norm.F|6107533|
|CASE_J03_trim_norm.R|6107533|
|CASE_J09_trim_norm.F|6999208|
|CASE_J09_trim_norm.R|6999208|
|CASE_J12_trim_norm.F|6058340|
|CASE_J12_trim_norm.R|6058340|
|CASE_J13_trim_norm.F|7623465|
|CASE_J13_trim_norm.R|7623465|
|CON_J02_trim_norm.F|6292686|
|CON_J02_trim_norm.R|6292686|
|CON_J05_trim_norm.F|6376901|
|CON_J05_trim_norm.R|6376901|
|CON_J10_trim_norm.F|7106205|
|CON_J10_trim_norm.R|7106205|
|SE_J01_trim_norm.F|6810374|
|SE_J01_trim_norm.R|6810374|
|SE_J04_trim_norm.F|6389705|
|SE_J04_trim_norm.R|6389705|
|SE_J07_trim_norm.F|6165907|
|SE_J07_trim_norm.R|6165907|

```bash
fastqc -f fastq *.fq
multiqc -o multiqc .
```

```bash
# use conda to activate pear env
# conda activate pear
# set base variable as merge directory for output
base="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/merge/"
# declare a string array with all of the filename prefixes
declare -a StringArray=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
# for loop to use all file prefixes and run pear aligner
# -j: threads
# -f: forward reads
# -r: reverse reads
# -o: output directory and filename
for i in "${StringArray[@]}"
do
pear -j 20 -f "$i"_trim.F.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq -r "$i"_trim.R.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq -o "$base""$i"_trim_norm_merge.fq.gz
done
```

Interesting note that pear automatically adds a "assembled.fastq" to the end of your putput files. I did also clear up these directories and gzipped them to help save space

Time to count all of the lines within the merge fastq files

```bash
# foor loop that
for fq in *assembled.fastq
do
echo $fq
cat $fq | echo $((`wc -l`/4))
done
```

|Sample Name|Read Count|
|-----------|----------|
|CA_J06_merge|4552928|
|CA_J08_merge|5001324|
|CA_J11_merge|5475886|
|CA_J18_merge|6678405|
|CASE_J03_merge|5402832|
|CASE_J09_merge|6278656|
|CASE_J12_merge|5328636|
|CASE_J13_merge|6816319|
|CON_J02_merge|5753317|
|CON_J05_merge|5572087|
|CON_J10_merge|6397342|
|SE_J01_merge|5977516|
|SE_J04_merge|5771304|
|SE_J07_merge|5341931|

### Cluster

First we need to bring all the normalized merged filed together

```bash
cat *assembled.fastq >> all.norm.merge.assembled.fastq
```

Next we need to cluster all the reads with cd-hit-est.
```bash
# cd hit is a fairly simple program that will cluster reads reducing the amount of information we are processing as this "seed"
# -i: input file
# -o: output file name
# -c: sequence identity threshold
# -n: word length
# -d: length of description in .clstr file, default 20 if set to 0, it takes the fasta defline and stops at first space
# -M: memory limit in MB
# -T: threads
# use conda to create cdhit env
# mamba create -n cdhit cd-hit
cd-hit-est -i all.norm.merge.assembled.fastq -o cDNA.c80.cluster -c 0.80 -n 5 -d 0 -M 26000 -T 10
```

Convert fastq to fasta file for use in itero

```bash
# seqtk is a program that has a ton of great tools to manipulate fastq and fasta files
# -N: remove N'sbow
# # use conda to create seqkit env
# mamba create -n seqkit seqkit
seqtk seq -N all.norm.merge.assembled.fastq > all.norm.merge.assembled.fasta
```

```bash
cd-hit-est -i all.norm.merge.assembled.fastq -o ~/repos/BIO594_work/course_project/data/cDNA/cluster/cDNA.norm.assembled.c95.cluster -c 0.80 -n 10 -d 0 -M 16000 -T 10
```

## Itero assembly

Create input config file for the program to use which has the following structure

```bash
[reference]
/home/jgreen/repos/BIO594_work/course_project/data/cDNA/merge_norm/all.norm.merge.assembled.fa.gz

[individuals]
CA_J06:/home/jgreen/repos/BIO594_work/course_project/data/gDNA/norm/CA_J06/
```

```bash
#mamba create -n itero itero
itero assemble local --config case.conf --output local --local-cores 10 --iterations 6
```

This has failed alot. There seems to be an issue with bam headers and potentially having to many which take up a ton of space. So I went looking for a memory input variable as there isn't a flag in the program options. The only one I could find for the program itero is in its itero.conf file within the miniconda envs itero folder. I changed the memory under spades from 2 to 50.

```bash
bedtools:$CONDA_PREFIX/bin/bedtools
bwa:$CONDA_PREFIX/bin/bwa
gawk:$CONDA_PREFIX/bin/gawk
grep:$CONDA_PREFIX/bin/grep
samtools:/usr/local/bin/samtools
spades:$CONDA_PREFIX/bin/spades.py

[spades]
kmer:33
coverage_cutoff:5
# 'memory' is RAM used per proc in GB
memory:50
```

## Trinity

To run trinity, I'll cat the forward and reverse files into one file respectively and use then to build a de novo transcriptome.

```bash
cat *.F*.fq > cDNA.F.fq
cat *.R*.fq > cDNA.R.fq
```

Run trinity

```bash
Trinity --no_normalize_reads --seqType fq --max_memory 100G --left cDNA.F.fq.gz  --right cDNA.R.fq.gz --CPU 10
```

Once trinity assembly is finished run the Trinitysta.pl

```bash
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  348348
Total trinity transcripts:      648435
Percent GC: 38.76

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 5280
        Contig N20: 3537
        Contig N30: 2552
        Contig N40: 1852
        Contig N50: 1312

        Median contig length: 398
        Average contig: 765.37
        Total assembled bases: 496294732


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 4378
        Contig N20: 2669
        Contig N30: 1654
        Contig N40: 1014
        Contig N50: 675

        Median contig length: 329
        Average contig: 556.54
        Total assembled bases: 193868658
```

Run cd-hit-est on this assembly


```bash
cd-hit-est -i trinity_out_dir.Trinity.fasta -o trinity.transcripts.cluster90.fasta -n 10 -d 0 -M 16000 -T 10
```

Once cluserting the trinity assembly is finished run the Trinitysta.pl

```bash
TrinityStats.pl trinity.transcripts.cluster90.fasta

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  323613
Total trinity transcripts:      434806
Percent GC: 38.29

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 4719
        Contig N20: 3040
        Contig N30: 2043
        Contig N40: 1359
        Contig N50: 894

        Median contig length: 358
        Average contig: 638.23
        Total assembled bases: 277504415


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 4462
        Contig N20: 2768
        Contig N30: 1748
        Contig N40: 1091
        Contig N50: 721

        Median contig length: 336
        Average contig: 575.23
        Total assembled bases: 186151685
```

### Assessing Trinity Assembly

#### Detonate

```bash
rsem-eval-estimate-transcript-length-distribution input.fasta parameter_file
rsem-eval-calculate-score -p 5 --transcript-length-parameters trinity.rsem.eval.ETLD.param --paired-end cDNA.F.fq.gz cDNA.R.fq.gz trinity.transcripts.cluster90.fasta trinity_cluster 638
```
Output from Rsem-eval score file

```bash
Score	-17910641055.22
BIC_penalty	-3982401.69
Prior_score_on_contig_lengths_(f_function_canceled)	-447798.78
Prior_score_on_contig_sequences	-384702805.70
Data_likelihood_in_log_space_without_correction	-17523466405.51
Correction_term_(f_function_canceled)	-1958356.47
Number_of_contigs	434806
Expected_number_of_aligned_reads_given_the_data	45800154.04
Number_of_contigs_smaller_than_expected_read/fragment_length	335638
Number_of_contigs_with_no_read_aligned_to	43536
Maximum_data_likelihood_in_log_space	-17498868545.99
Number_of_alignable_reads	52009248
Number_of_alignments_in_total	99279348
```

Output from RSEM-eval cluster gene results

```bash
gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
TRINITY_DN0_c0_g1_i3	TRINITY_DN0_c0_g1_i3	974.00	763.32	528.54	15.84	15.10
TRINITY_DN0_c0_g1_i5	TRINITY_DN0_c0_g1_i5	934.00	723.32	554.87	17.55	16.73
TRINITY_DN0_c0_g1_i9	TRINITY_DN0_c0_g1_i9	1011.00	800.32	122.59	3.50	3.34
TRINITY_DN0_c1_g1_i10	TRINITY_DN0_c1_g1_i10	4244.00	4033.32	3880.42	22.01	20.99
TRINITY_DN0_c2_g1_i2	TRINITY_DN0_c2_g1_i2	498.00	287.68	10.31	0.82	0.78
TRINITY_DN0_c2_g1_i4	TRINITY_DN0_c2_g1_i4	389.00	179.91	10.00	1.27	1.21
TRINITY_DN0_c2_g1_i7	TRINITY_DN0_c2_g1_i7	241.00	46.28	0.00	0.00	0.00
TRINITY_DN0_c2_g1_i8	TRINITY_DN0_c2_g1_i8	588.00	377.47	87.68	5.31	5.07
```

#### Transrate

```bash
transrate --assembly transcripts.fa --left left.norm.fq --right right.norm.fq --threads 10 --output transrate_oases 
```
|assembly|n_seqs|smallest|largest|n_bases|mean_len|n_under_200|n_over_1k|n_over_10k|n_with_orf|mean_orf_percent|n90|n70|n50|n30|n10|gc|bases_n|proportion_n|fragments|fragments_mapped|p_fragments_mapped|good_mappings|p_good_mapping|bad_mappings|potential_bridges|bases_uncovered|p_bases_uncovered|contigs_uncovbase|p_contigs_uncovbase|contigs_uncovered|p_contigs_uncovered|contigs_lowcovered|p_contigs_lowcovered|contigs_segmented|p_contigs_segmented|score|optimal_score|cutoff|weighted|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|/home/jgreen/repos/BIO594_work/course_project/data/assembly/trinity/trinity_out_dir.Trinity.fasta|648435|177|40397|496294732|765.32192|173|120054|768|100909|48.46403|299|595|1313|2553|5281|0.38756|0|0.0|90243319|32124228|0.35597|8320550|0.0922|23803678|0|186218129|0.37522|470161|0.72507|648435|1.0|648435|1.0|68753|0.10603|0.00958|0.02607|0.13427|0.0|

#### RNAquast

```bash
rnaQUAST.py -c trinity_out_dir.Trinity.fasta -1 cDNA.F.fq.gz -2 cDNA.R.fq.gz --busco ~/databases/mollusca_odb10/ -t 10 -o trinity_rnaquast
```

SHORT SUMMARY REPORT 

METRICS/TRANSCRIPTS                                    trinity_out_dir.Trinity  

 == BASIC TRANSCRIPTS METRICS == 
Transcripts                                            648435                   
Transcripts > 500 bp                                   249381                   
Transcripts > 1000 bp                                  120165

 == BUSCO METRICS == 
Complete                                               97.658                   
Partial                                                1.228

### gDNA reads

Symbolic link for reads in the gDNA folder for gDNA reads

```bash
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture2-2-Capture3-1_R1_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture2-2-Capture3-1_R2_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture1-1-Capture2-1_R1_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture1-1-Capture2-1_R2_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture1-2-Capture3-2-Capture4-1_R1_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture1-2-Capture3-2-Capture4-1_R2_001.fastq.gz 
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture2-3-Capture3-3-Capture4-2_R1_001.fastq.gz
ln -s /RAID_STORAGE2/Raw_Data/CASE/Block11/Capture2-3-Capture3-3-Capture4-2_R2_001.fastq.gz
```

Run fastqc for all gDNA reads that are not de multiplexed. We will be keeping these reads demultiplexed for assembly

```bash
fastqc -f fastq *.fastq
multiqc -o multiqc .
```

Next we are going to trim each sample with fastp

```bash
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/CA_J06.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/CA_J06.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/Capture1-1-Capture2-1_R1_001.fastq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/Capture1-1-Capture2-1_R1_001.fastq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/CA_J06_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/CA_J06_trim.R.fq.gz
```

Run fastqc on all trimmed samples

```bash
fastqc -f fastq *.fastq
multiqc -o multiqc .
```

These reads were already demultiplexed so I symbolic linked them to a demux directory

```bash
ln -s /RAID_STORAGE2/Shared_Data/CASE/Block11/capture_1-1_2-1/
ln -s /RAID_STORAGE2/Shared_Data/CASE/Block11/capture_1-1_2-1/
ln -s /RAID_STORAGE2/Shared_Data/CASE/Block11/capture_1-1_2-1/
ln -s /RAID_STORAGE2/Shared_Data/CASE/Block11/capture_1-1_2-1/
```

Time to make a for loop through a list of the names of the demultiplexed file types for each capture

```bash
declare -a StringArray2=("CA_J06" "CASE_J03" "CON_J02" "CON_J05" "SE_J01" "SE_J04")

for i in "${StringArray2[@]}"
do
echo $i
done

for i in "${StringArray2[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-1_2-1/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-1_2-1/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.R.fq.gz
done
```

Lets make a fastp.sh script and set our stringarray for every different capture combination

```bash
declare -a StringArray3=("CA_J18" "CASE_J12" "CASE_J13" "CON_J15" "IS_01" "SE_J14")

for i in "${StringArray3[@]}"
do
echo $i
done

for i in "${StringArray3[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-2_3-2/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-2_3-2/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.R.fq.gz
done

declare -a StringArray4=("CA_J08" "CA_J11" "CASE_J09" "CON_J10" "SE_J07")

for i in "${StringArray4[@]}"
do
echo $i
done

for i in "${StringArray4[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_2-2_3-1/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_2-2_3-1/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.R.fq.gz
done

declare -a StringArray5=("CA_J18G" "CASE_J13G" "IS_02" "IS_03" "IS_04" "SE_J14G")

for i in "${StringArray5[@]}"
do
echo $i
done

for i in "${StringArray5[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_2-3_3-3/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_2-3_3-3/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.R.fq.gz
done
```

Time to run fastqc on all the different trimmed FQ files

```bash
fastqc -f fastq *.fastq
multiqc -o multiqc .
```

```bash
for fq in *.fq.gz
do
echo $fq
cat $fq | echo $((`wc -l`/4))
done
```

|Sample Name|Read Count|
|-----------|----------|
|CA_J06_trim.F.fq.gz|7771638|
|CA_J06_trim.R.fq.gz|7771638|
|CA_J08_trim.F.fq.gz|7195617|
|CA_J08_trim.R.fq.gz|7195617|
|CA_J11_trim.F.fq.gz|10479050|
|CA_J11_trim.R.fq.gz|10479050|
|CA_J18G_trim.F.fq.gz|9018258|
|CA_J18G_trim.R.fq.gz|9018258|
|CA_J18_trim.F.fq.gz|11015798|
|CA_J18_trim.R.fq.gz|11015798|
|CASE_J03_trim.F.fq.gz|7706833|
|CASE_J09_trim.F.fq.gz|8245369|
|CASE_J09_trim.R.fq.gz|8245369|
|CASE_J12_trim.F.fq.gz|9075471|
|CASE_J12_trim.R.fq.gz|9075471|
|CASE_J13G_trim.F.fq.gz|8411467|
|CASE_J13G_trim.R.fq.gz|8411467|
|CASE_J13_trim.F.fq.gz|7845068|
|CASE_J13_trim.R.fq.gz|7845068|
|CON_J02_trim.F.fq.gz|9573138|
|CON_J02_trim.R.fq.gz|9573138|
|CON_J05_trim.F.fq.gz|5699313|
|CON_J05_trim.R.fq.gz|5699313|
|CON_J10_trim.F.fq.gz|7213905|
|CON_J10_trim.R.fq.gz|7213905|
|CON_J15_trim.F.fq.gz|7299912|
|CON_J15_trim.R.fq.gz|7299912|
|IS_01_trim.F.fq.gz|11652092|
|IS_01_trim.R.fq.gz|11652092|
|IS_02_trim.F.fq.gz|5855570|
|IS_02_trim.R.fq.gz|5855570|
|IS_03_trim.F.fq.gz|8604171|
|IS_03_trim.R.fq.gz|8604171|
|IS_04_trim.F.fq.gz|8984037|
|IS_04_trim.R.fq.gz|8984037|
|SE_J01_trim.F.fq.gz|8052394|
|SE_J01_trim.R.fq.gz|8052394|
|SE_J04_trim.F.fq.gz|5840686|
|SE_J04_trim.R.fq.gz|5840686|
|SE_J07_trim.F.fq.gz|6902310|
|SE_J07_trim.R.fq.gz|6902310|
|SE_J14G_trim.F.fq.gz|8850953|
|SE_J14G_trim.R.fq.gz|8850953|
|SE_J14_trim.F.fq.gz|8219419|
|SE_J14_trim.R.fq.gz|8219419|

```bash
baseg="/home/jgreen/repos/BIO594_work/course_project/data/gDNA/norm/"
#readsg="/home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/"
declare -a StringArray=("CA_J06_trim" "CA_J08_trim" "CA_J11_trim" "CA_J18G_trim" "CA_J18_trim" "CASE_J03_trim" "CASE_J09_trim" "CASE_J12_trim" "CASE_J13G_trim" "CASE_J13_trim" "CON_J02_trim" "CON_J05_trim" "CON_J10_trim" "CON_J15_trim" "IS_01_trim" "IS_02_trim" "IS_03_trim" "IS_04_trim" "SE_J01_trim" "SE_J04_trim" "SE_J07_trim" "SE_J14G_trim" "SE_J14_trim")
for i in "${StringArray[@]}"
do
insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 30 --left "$i".F.fq.gz --right "$i".R.fq.gz --pairs_together --PARALLEL_STATS --CPU 10 --output "$baseg"
done
```

```bash
baseg="/home/jgreen/repos/BIO594_work/course_project/data/gDNA/norm/"
declare -a StringArray=("CA_J06" "CA_J08" "CA_J11" "CA_J18G" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13G" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "CON_J15" "IS_01" "IS_02" "IS_03" "IS_04" "SE_J01" "SE_J04" "SE_J07" "SE_J14G" "SE_J14")
for i in "${StringArray[@]}"
do
	cd $i
	for fq in *READ*.fq
	do
		echo $fq
		cat $fq | echo $((`wc -l`/4))
	done
	cd $baseg
done
```

|Sample Name|Read Count|
|-----------|----------|
|CA_J06-READ1.fq|6474353|
|CA_J06-READ2.fq|6474353|
|CA_J08-READ1.fq|6102510|
|CA_J08-READ2.fq|6102510|
|CA_J11-READ1.fq|8737519|
|CA_J11-READ2.fq|8737519|
|CA_J18G-READ1.fq|7267360|
|CA_J18G-READ2.fq|7267360|
|CA_J18-READ1.fq|9085969|
|CA_J18-READ2.fq|9085969|
|CASE_J03-READ1.fq|6314595|
|CASE_J03-READ2.fq|6314595|
|CASE_J08-READ1.fq|6865755|
|CASE_J08-READ2.fq|6865755|
|CASE_J12-READ1.fq|7642087|
|CASE_J12-READ2.fq|7642087|
|CASE_J13G-READ1.fq|6845444|
|CASE_J13G-READ2.fq|6845444|
|CASE_J13-READ1.fq|6647585|
|CASE_J13-READ2.fq|6647585|
|CON_J02-READ1.fq|7728305|
|CON_J02-READ2.fq|7728305|
|CON_J05-READ1.fq|4723536|
|CON_J05-READ2.fq|4723536|
|CON_J10-READ1.fq|6203104|
|CON_J10-READ2.fq|6203104|
|CON_J15-READ1.fq|6237440|
|CON_J15-READ2.fq|6237440|
|IS_01-READ1.fq|9307311|
|IS_01-READ2.fq|9307311|
|IS_02-READ1.fq|4544553|
|IS_02-READ2.fq|4544553|
|IS_03-READ1.fq|6594197|
|IS_03-READ2.fq|6594197|
|IS_04-READ1.fq|6888736|
|IS_04-READ2.fq|6888736|
|SE_J01-READ1.fq|6618882|
|SE_J01-READ2.fq|6618882|
|SE_J04-READ1.fq|4935837|
|SE_J04-READ2.fq|4935837|
|SE_J07-READ1.fq|5755544|
|SE_J07-READ2.fq|5755544|

## Mapping gDNA reads to de novo transcriptome

```bash
bwa index -a bwtsw trinity.transcripts.cluster90.fasta trinity.bwa
bwa mem -t 20 trinity.transcripts.cluster90.fasta gDNA.F.fq gDNA.R.fq | samtools view -bu - | samtools sort -@4 - -o trinity.c90.sorted.bam
```

## Trinity transcriptome guided assembly

```bash
 Trinity --genome_guided_bam trinity.c90.sorted.bam \
         --genome_guided_max_intron 5000 \
         --max_memory 50G --CPU 20
```
Unfortunately the project will have to end abruptly here. This specific line of code generated ~ 1 TB of data overnight. We are reconsidering how to subset data and creating upperthresholds and standard practice for analyzing this many samples.

```bash
~/miniconda3/envs/trinity/util/misc/process_GMAP_alignments_gff3_chimeras_ok.pl \
     --genome trinity_out_dir.Trinity.fasta \
     --transcripts trinity_out_dir/<filename>.fasta \
     --SAM | samtools view -Sb | samtools sort -o trinity-GG.gmap.bam
```