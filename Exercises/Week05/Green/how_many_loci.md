# Take-home Exercise for Week 5

## How many loci are there?

I've placed some files in the `/home/BIO594/DATA/Week5/` directory.  They are ddRADseq files that have already been demultiplexed.  

Using your new RADSeq assembly skills:

* Post a markdown document to this directory showing your efforts to answer the question: 
	* How many loci are there in the data set?

* Post a `reference.fasta` file with your assembled reference

* You can use the material from class or anything else you might find at [dDocent.com](dDocent.com)

## Process of finding how many loci

Create set of uniq reads with counts for each individual

```{bash eval=FALSE}
 ls *.F.fq.gz > namelist
 sed -i'' -e 's/.F.fq.gz//g' namelist
 AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
 AWK2='!/>/'
 AWK3='!/NNN/'
 PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
 
 cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
 cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
 cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
 ```

Sum up the number the within individual converage level of unique reads in our dataset

```{bash eval=FALSE}
cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
```

Use mawk to query the first column and select data above a certain copy number (from 2-20) and prints that to a file.

```{bash eval=FALSE}
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile
```

Look at the contents of uniqseq.data

```{bash eval=FALSE}
more uniqseq.data
```

Plot uniseq.data

```{bash eval=FALSE}
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
     130000 +-------------------------------------------------------------------------------------------------------+   
            |           +          +           +          +           +          +           +          +           |   
     125000 |**                                                                                                   +-|   
            |  ******                                                                                               |   
            |        ******                                                                                         |   
     120000 |-+            *****                                                                                  +-|   
            |                   ******                                                                              |   
     115000 |-+                       ******                                                                      +-|   
            |                               ******                                                                  |   
     110000 |-+                                   *****                                                           +-|   
            |                                          ******                                                       |   
     105000 |-+                                              ******                                               +-|   
            |                                                      ******                                           |   
            |                                                            ******                                     |   
     100000 |-+                                                                ****                               +-|   
            |                                                                      ***                              |   
      95000 |-+                                                                       ****                        +-|   
            |                                                                             ******                    |   
      90000 |-+                                                                                 *****             +-|   
            |                                                                                        **             |   
            |                                                                                          ****         |   
      85000 |-+                                                                                            *****  +-|   
            |           +          +           +          +           +          +           +          +       *** |   
      80000 +-------------------------------------------------------------------------------------------------------+   
            2           4          6           8          10          12         14          16         18          20  
                                                            Coverage     


Using cutoff value 4

```{bash eval=FALSE}
parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv
```

Data reduced to 12699 sequences with cutoff value 4

Further reduce data by restricting the number of different individuals a seuence appears

```{bash eval=FALSE}

for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done
 
cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data

rm ufile
```

plot the values from individual restriction

```{bash eval=FALSE}

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
set xrange [2:20]
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF

```

Graph of Reduction

                                                                                                                       
                                                                                                                        
                       Number of Unique Sequences with More than X Coverage (Counted within individuals)                
     11000 +--------------------------------------------------------------------------------------------------------+   
           |           +          +           +           +          +           +           +          +           |   
           |                                                                                                        |   
           |                                                                                                        |   
     10000 |*+                                                                                                    +-|   
           | *                                                                                                      |   
           |  **                                                                                                    |   
           |    *                                                                                                   |   
      9000 |-+   *                                                                                                +-|   
           |      **                                                                                                |   
           |        *                                                                                               |   
      8000 |-+       **                                                                                           +-|   
           |           **                                                                                           |   
           |             ***                                                                                        |   
           |                ***                                                                                     |   
      7000 |-+                 **                                                                                 +-|   
           |                     ***                                                                                |   
           |                        ***                                                                             |   
           |                           ***                                                                          |   
      6000 |-+                            ***                                                                     +-|   
           |                                 ****                                                                   |   
           |                                     *****                                                              |   
           |           +          +           +       *** +          +           +           +          +           |   
      5000 +--------------------------------------------------------------------------------------------------------+   
           2           4          6           8           10         12          14          16         18          20  
                                                                       Coverage  

Using cutoff value 4

```{bash eval=FALSE}
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs
```
We have reduced the number of sequences to 7989 in uniq.k.4.c.4.seqs

Convert these sequences into fasta files

```{bash eval=FALSE}
cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
```

Extract forward reads

```{bash eval=FALSE}
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
```
Cluster forward reads at 88% similarity with cd hit est.

```{bash eval=FALSE}
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1
```

Convert CD hit output into format for furst step of rainbow

```{bash eval=FALSE}
mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
```
Count the number of clusters

```{bash eval=FALSE}
cut -f2 rcluster | uniq | wc -l
```
We have approximately 1234 clusters

Now we split these clusters in smaller clusters representing significant variants and count lines

```{bash eval=FALSE}
rainbow div -i rcluster -o rbdiv.out
wc -l rbdiv.out
```
Output shows 7989 significant variants

Lets alter some parameters for frequency (-f) and the minimum number of allelles split (-K) and count the difference in significant variants

```{bash eval=FALSE}
rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
wc -l rbdiv.out
```
Output shows 7989 significant variants, so this had no effect. Lets try more strict parameters to ellicit an effect

```{bash eval=FALSE}
rainbow div -i rcluster -o rbdiv.out -f 0.8 -K 5
wc -l rbdiv.out
```

No effect either way

Using paried end reads to merge clusters with a cutoff value of 2

```{bash eval=FALSE}
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
```

This next section looks at all the contigs assembled for a cluster.  If any of the contigs contain forward and PE reads, then that contig is output as optimal. If no overlap contigs exists (the usual for most RAD data sets), then the contig with the most assembled reads PE (most common) is output with the forward read contig with a 10 N spacer.   If two contigs have equal number of reads, the longer contig is output. 

I have questions for this section

```{bash eval=FALSE}
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

Checking for overlap (did not do this but want to come back aroudn to this section)

```{bash eval=FALSE}
seqtk seq -r rainbow.fasta > rainbow.RC.fasta
mv rainbow.RC.fasta rainbow.fasta

#The rainbow assembly is checked for overlap between newly assembled Forward and Reverse reads using the software PEAR
sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f1 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.F.fq
sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f2 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.R.fq

seqtk seq -r ref.R.fq > ref.RC.fq
mv ref.RC.fq ref.R.fq
LENGTH=$(mawk '!/>/' rainbow.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(( $LENGTH * 5 / 4))
	
pearRM -f ref.F.fq -r ref.R.fq -o overlap -p 0.001 -j $NUMProc -n $LENGTH

rm ref.F.fq ref.R.fq

mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.assembled.fastq > overlap.fasta
mawk '/>/' overlap.fasta > overlap.loci.names
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.forward.fastq > other.F
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.reverse.fastq > other.R
paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNN/g' > other.FR

cat other.FR overlap.fasta > totalover.fasta

rm *.F *.R

fi
if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" ]]; then
	cp uniq.fasta totalover.fasta
fi
cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T 0 -c $simC

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta
```

Align and cluster sequences

```{bash eval=FALSE}
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9
```

Remake the reference with a cutoff of 20 copies of a unique sequence to use for assembly and a final clustering value of 90%.

```{bash eval=FALSE}
bash remake_reference.sh 4 4 0.90 PE 2
```

This script uses different loops to assemble references from an interval of cutoff values and c values from 0.8-0.98.

```{bash eval=FALSE}
bash ReferenceOpt.sh 4 8 4 8 PE 16 &
disown -a
```