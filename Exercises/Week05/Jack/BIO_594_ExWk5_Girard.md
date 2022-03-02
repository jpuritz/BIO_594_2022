# Week 5 Homework Assignment RADseq alignment, reference, and Loci
Jack Girard | 01-Mar-2022

## Final Output

Loci #1234

### Genomics parameters

Coverage: 8

Individuals: 4

Similarity Percentages Tested: 80 & 90

## Resouces and Goals

[Goals Link](https://github.com/jpuritz/BIO_594_2022/blob/main/Exercises/Week05/Exercise.md)

[Class Notes Link](https://docs.google.com/document/d/1s3NWW4IjVGKv_rDxtDpzRr3qBSJEbq6ThfSOF0s9Gek/edit#)

Data located in `/home/BIO594/DATA/Week5/`

* Looks like we have 5 populations A–D each with 20 individuals sampled in each.

* No barcodes file

* Data has already been demultiplexed

## Initial environment

>~~Creating a symbolic link to the data~~ Was unable to properly manipulate and reference the files. The workflow is easier if copied.

Copying data to a directory called loci

Starting directory on kitt `/home/jgirard/wk5`

````
```bash
mkdir goal
cd goal

cp -r  /home/BIO594/DATA/Week5/ ./loci 

ls loci # files present
```
````
Inspecting files 

````
```bash
zcat loci/PopA_001.F.fq.gz | head -15
zcat loci/PopA_002.F.fq.gz | head -15
zcat loci/PopD_002.F.fq.gz | head -15
zcat loci/PopB_005.F.fq.gz | head -15
```
````
sequenes appear to be nearly identical. The all start with AATT and end with GCC. Some quick googling shows these to be restriction sites 

We will be using a conda environment with `ddocent` loaded and the `process_radtags` function.

````
```bash
conda info --env #Ref.Ex envi has ddocent loaded
conda activate Ref.Ex
```
````

## Assembling a Reference with the dDocent pipeline

### Creating a list of unique reads

Creating a name list | setting variables for the environment | creating files with those variables

````
```bash

cd loci

ls *.F.fq.gz > namelist #creating name list
sed -i'' -e 's/.F.fq.gz//g' namelist # find and replace within namelist

#below this point setting variables 

AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

#below this point creating new files using the variables defined above

cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward" #creates foreward reads removing quality scores
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse" #creates PE reads same way as above
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs" #puts both reads together with 10 N's between them, finds the unique reads, counts them
 ```
 ````

Each sample now has a `.forward`, `.reverse`, and `.uniq.seqs` file

>Due to the nature of RADseq usually reads with _low copy_ numbers can be eliminated. This is because sequences with _low coverage_ tend to be sequence errors

 Summing up coverage within individuals for unique seqs.

BASH for loop from Ref.Ex tutorial

````
```bash

cat *.uniq.seqs > uniq.seqs

for i in {2..20};
do
echo $i >> pfile
done

cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data

$rm pfile
```
````

`uniqseq.data`, `uniq.seqs` files created.

This data can be plotted directly to the terminal to view the coverage of Unique Sequenes

````
```bash
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
````

![](2022-03-01-21-22-28.png)

Interpretting the graph it looks like thresholding uniq seqs with coverage of 8 will preseve most most of sequences while limiting risk.

**Choosing a coverage cutoff value:** preserve maximum diversity retain minimal sequence errors.

````
```bash
parallel --no-notice -j 8 mawk -v x=8 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv

wc -l uniqCperindv
```
````
>12442 sequences at greater than 8 coverage

The x = after the `mawk -v` sets the cutoff coverage value

Subsetting sequences based on the number of different individuals in which a particular squence occures. We are choosing 2 to see the range of data 

````
```bash
for ((i = 2; i <= 2; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data

rm ufile
```
````
`i <= 2` is were the number of individualsis specified.

Replotting data as the number of sequences occuring within individuals at greater than 8 coverage

````
```bash
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
````

![](2022-03-01-22-38-01.png)

Appears to be a rapid decline as individuals increases. Convervatively 4 individuals may be the best bet without loosing more than half of the sequences.

Rerunning individuals subset for 4 individuals.

````
```bash
for ((i = 2; i <= 4; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data

rm ufile
```
````

using `mawk` to specify a coverage cut off of 4 – this time for the minimum number of individuals

````
```bash
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.8.seqs #Note the change from "<" to ">"

wc -l uniq.k.4.c.8.seqs #checking samples left
```
````
>7743 seqs left

### converting seqs back to `fasta` format

With data refind reformatting

````
```bash
cut -f2 uniq.k.4.c.8.seqs > totaluniqseq  #subset the 2nd col file and create new file
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
```
````
### Assembling Refernce Contigs

First extracting the forward reads. Remember the foreward and reverse reads were seperated by 10 N's

````
```bash
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
```
````
Attempting to break this down

* -e — `man sed` says this tells `sed` it is a script "add script to commands to be executed"
* \' \' — signifies string
* s — search?
* /NNNNNNNNNN/ — Find this pattern
* \t — replace with a `TAB`
* /g — perform globally, change all occurances
* input file uniq.fasta
* | — pipe to `cut` 
* -f1 — take the first column
* output as uniq.F.fasta 

#### Clustering and Alignment

>Best to do this starting with "low" similarity (e.g 80%) and let other functions futher subset the data later.

Using `CD-HIT` 

Lab [Link](http://weizhong-lab.ucsd.edu/cd-hit/)

[Github](https://github.com/weizhongli/cdhit)

[User .pdf](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.pdf)

````
```bash
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1
```
````

Breaking it down
> To see `cd-hit` documentation type in the command with no options set.

* cd-hit-est — the "est" is specific for neucleotides. Otherwise  the progam assumes it is working with peptides
* -i — flag for input file
* -o — flag for output file
* -c — flag for proportion of similarity 1.00 = 100%
* -T — Sets the number of threads setting the value to 0 uses all available.
* -M — Memory, specifies the amount of RAM to use. Setting to 0 uses all available.
* -g — can be 1 or 0 (default). Changes how the algorithm clusters. 1 is slower yet more accurate according to the function documentation

Next take `cd-hit` output and convert it into a format for `rainbow` to handle.

````
```bash
mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids

paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq

sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
```
````
Inspect files with `more` `head` and `tail`. To determine the final number of clusters `cut` `uniq` and `wc` can be piped together

````
```bash
cut -f2 rcluster | uniq | wc -l #col 2 has the cluster ID
```
````
>1234 clusters RAD loci found

##### Clustering by alleles

With the RAD loci found now we can begin to cluster further (i.e. alleles)

````
```bash

rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
```
````
Breaking it down

* -f — minimum frequeny needed to divided an allele into a unique cluster
* -K – controls minimum number of alleles to split regarless of frequency

##### Merging Divided Clusters

>This partially serves as a QA/QC check. Theory suggests that Paired End reads should be very similar (i.e. the reading directions should not matter) Those with very different seqs are more likely to be "repetitive regions"

````
```bash

rainbow merge -o rbasm.out -a -i rbdiv.out -r 2

```
````
* -r — minimum number of reads to assemble \[5\]

#### Optimal & Suboptimal Conitigs

dDocent will do this however, Ref.Ex skips this step. Below is the code and description directly from the Ref.Ex activity.

>Now, this looks a bit complicated, but it's performing a fairly simple algorithm.
 First, the script looks at all the contigs assembled for a cluster.  If any of the contigs contain forward and PE reads, then that contig is output as optimal.     
 If no overlap contigs exists (the usual for most RAD data sets), then the contig with the most assembled reads PE (most common) is output with the forward read contig with a 10 N spacer.
 If two contigs have equal number of reads, the longer contig is output.

 At this point, dDocent (version 2.0 and higher) will check for substantial overlap between F and PE reads in the contigs.

````
``bash
cat rbasm.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
if (NR == 1) e=$2;
else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/C/) clus=$2;
else if ($1 ~/L/) len=$2;
else if ($1 ~/S/) seq=$2;
else if ($1 ~/N/) freq=$2;
else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
}' > rainbow.fasta
```
````
With optimized contigs we will now increase the similarity threashold

````
```bash
cd-hit-est -i rainbow.fasta -o girard.reference.fasta -M 0 -T 0 -c 0.9
```
````
The number of clusters (loci) is reported after the program runs
>1234 clusters (loci)
No change with increased similarity threshold.

If others want to automate this process to check my data below is a shell script that will allow for the manipulation of the coverage, presence in individuals, and similarity percentage parameters

````
```bash
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh

bash remake_reference.sh 4 4 0.90 PE 2
```
````

## PyRAD Environment Set Up

Loading conda environment and creating pyrad parameters file

````
```bash
conda info --envs # checking environments. "pyrad" environment present

conda activate pyrad

pyrad -n # creates a file params.txt
```
````
Editing `params.txt` to include: 
* #6 restriction sites AATT, CCG
* #7 Setting parallel processors to 12
* #11 settng to paird ddRAD
* #12 set trim overhang --> preventing a pyrad specific bug in pyRAD 
* #13 setting to 20 individuals

