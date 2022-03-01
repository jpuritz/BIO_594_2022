# Take-home Exercise for Week 5

### How many loci are there?

##### TLDR: There are 1235 loci in this data. 

Using code from dDocent [site](http://www.ddocent.com/assembly/)

Make directory for this exercise 

```
mkdir Week5.Ex
cd Week5.Ex
```

We downloaded miniconda last week, so we can create and activate a dDocent conda environment to work in 

```
conda create -n ddocent_env ddocent

conda activate ddocent_env
```

Woohoo! Now we are working in the ddocent_env 

Copy the demultiplexed data into my own directory. The data is in the `/home/BIO594/DATA/Week5/` directory. 

```
cp /home/BIO594/DATA/Week5/*fq.gz .
less PopA_001.F.fq.gz # take a look at one of the files 
```

Now that we have the data, create a set of unique reads for each individual. First, set shell variables for awk and perl code

```
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
```

Create set of forward reads for each individual to sort fastq files and strip away quality scores. Then do same for reverse reads.

```
cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
```

Concatentates the forward & reverse reads, finds unique reads within that individual and counts occurences (coverage)

```
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
```

We want to get rid reads w/ low copy numbers as those are probably sequence errors and/or regions w/ very low coveraage. 

First, sum up number of within individual coverage level of unique reads

```
cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile
```

this code looks at the first column of uniq.seqs and selects the data above a certain copy number (in this case, 2-20) and prints to uniqseq.data file 

Check out the uniqseq.data file

```
more uniqseq.data 
2	126508
3	122492
4	121026
5	119357
6	117552
7	115633
8	113606
9	111449
10	109143
11	106710
12	104198
13	101560
14	98723
15	95844
16	92812
17	89829
18	86728
19	83535
20	80401
```

So I'm interpreting this file as 126,508 seqs have a copy number of 2? Can also plot it in terminal

```
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

```

Fun fun. Now we can choose a cutoff value that will capture as much diversity of the data as possible and eliminate seqs that are errors. 

Trying cutoff value of 4:

```
parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv
12699 uniqCperindv
```

We now have 12,699 sequences. We can further reduce by the number of individuals a sequence appears in

```
for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
rm ufile
```

Let's plot this one:

```
gnuplot << \EOF
> set terminal dumb size 120, 30
> set autoscale
> unset label
> set title "Number of Unique Sequences present in more than X Individuals"
> set xlabel "Number of Individuals"
> set ylabel "Number of Unique Sequences"
> plot 'uniqseq.peri.data' with lines notitle
> pause -1
> EOF

                                                                                                                        
                                                                                                                        
                                 Number of Unique Sequences present in more than X Individuals                          
     11000 +--------------------------------------------------------------------------------------------------------+   
           |            +            +            +             +            +            +            +            |   
           |                                                                                                        |   
           |*                                                                                                       |   
     10000 |-**                                                                                                   +-|   
           |   ***                                                                                                  |   
           |      ***                                                                                               |   
           |         **                                                                                             |   
      9000 |-+         ***                                                                                        +-|   
           |              ***                                                                                       |   
           |                 ****                                                                                   |   
      8000 |-+                   ***                                                                              +-|   
           |                        *****                                                                           |   
           |                             ******                                                                     |   
           |                                   *******                                                              |   
      7000 |-+                                        *******                                                     +-|   
           |                                                 *******                                                |   
           |                                                        ******                                          |   
           |                                                              *******                                   |   
      6000 |-+                                                                   ******                           +-|   
           |                                                                           **********                   |   
           |                                                                                     **********         |   
           |            +            +            +             +            +            +            +   ******   |   
      5000 +--------------------------------------------------------------------------------------------------------+   
           2            3            4            5             6            7            8            9            10  
                                                     Number of Individuals  
```

Let's cut off again at 4

```
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs
7989 uniq.k.4.c.4.seqs
```

Now we got 7,989 sequences. Cool so let's convert the seqs back to fasta format to move forward.

```
cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
zgrep -c ">" uniq.fasta
7989
```

At this point in the pipeline, dDocent usually checks for adapters, but since the data is simulated, we don't need to do that step.

Now we can assemble reference contigs. To begin, we can extract the forward reads from the fasta file 

```
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
```

This code replaces the 10N separator into tab character and then cuts it to split the files into forward reads. pretty smart

Now we will move on to full RAD assembly. Alignment clustering is done by CD-hit (used to be done by rainbow). 

```
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1

================================================================
Program: CD-HIT, V4.8.1 (+OpenMP), Feb 22 2022, 21:26:56
Command: cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g
         1

Started: Tue Mar  1 14:41:22 2022
================================================================
                            Output                              
----------------------------------------------------------------
total number of CPUs in the system is 80
Actual number of CPUs to be used: 80

total seq: 7989
longest and shortest : 116 and 116
Total letters: 926724
Sequences have been sorted

Approximated minimal memory consumption:
Sequence        : 1M
Buffer          : 80 X 12M = 972M
Table           : 2 X 16M = 33M
Miscellaneous   : 4M
Total           : 1013M

Table limit with the given memory limit:
Max number of representatives: 248606
Max number of word counting entries: 6715988

# comparing sequences from          0  to       7238
.......---------- new table with     1234 representatives
# comparing sequences from       7238  to       7989
....................---------- new table with        0 representatives

     7989  finished       1234  clusters

Approximated maximum memory consumption: 1014M
writing new database
writing clustering information
program completed !

Total CPU time 4.59
```

The code clusters the forward reads by 80% similarity. Convert the CD-hit output to match output of first phase of rainbow: 

```
mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
```

Check out the rcluster file

```
# output header
Read_ID	Cluster_ID	Forward_Read	Reverse_Read

more rcluster
1	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACATTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATAT
AGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
3269	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACGTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTACTCGTACCTGAACGTTTGTATGGATAT
AGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAAGC
4366	1	AATTCCAGGGCTGATATACCACCACACAAAAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCGAGATACTCCACTTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATAT
AGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
4385	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACATTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGCCTTTCGGTTACCGTAGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGACCTGAA
CTTTGTAGGATAAGGATCCACTCGGTACGGTCGCCTGCGCGAT

head rcluster 
1	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACATTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
3269	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACGTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTACTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAAGC
4366	1	AATTCCAGGGCTGATATACCACCACACAAAAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCGAGATACTCCACTTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
4385	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACATTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGCCTTTCGGTTACCGTAGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGACCTGAACTTTGTAGGATAAGGATCCACTCGGTACGGTCGCCTGCGCGAT
5924	1	AATTCCAGGGCTGATATACCACCACCCAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACGTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTACTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
6907	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACGTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
722	1	AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACGTTTATTTCATACGCTACTGTGCGACAGCCCTTT	CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTACTCGTACCTGAACGTTTGTATGGTTATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAAGC
1488	2	AATTCTTCGTTTGTCAAACGTTCAGGACGAAGGGTAACTGGAAGCATCGACGAGCTGAGATCTAGACTGTTCCAATTGAAGCTGAGCATCCACATTTTCTGCCGTAAGAACGCGCG	CGGTTCTGCCGCTGTGCGTGTGCGAGAGCGAGTACTCCACCCACCAAGCACAATCATACGTCCGCCTCGACGACTCCGCCCTATCCTTAGAAGGTGCGATGGGGTGCACTAGCGGGGCAAGCGCA
1708	2	AATTCTTCGTTTGTCAAACGTTCAGGACGAAGGGTAACTGGAAGCATCGACGAGCTGAGATCTAGACTGTTCCAATTGAAGCTGAGCATCCACATTTTCTGCCGTAAGAACGCGCG	CGGTTCTGCGCATGTGCGTGTGCGAGAGCGAGTACTCCACCCACCAAGCACAATCATACGTCCGCCTGAGCGACTCCGCCCTATCCTTAGAAGGTGCGATGGAGTGACTAGCGGGGCAAGCGCAT
2150	2	AATTCTTCGTTTGTCAAACGTTCAGGACGAAGGGTGACTGGAAGCATCGACGAGCTGAGATCTAGACTGTTCCAATTGAAGCTGAGCATCCACATTTTCTGCCGTAAGAACGCGCG	CGGTTCTGCCGCTGTGCGTGTGCGAGAGCGAGTACTCCACCCACCAAGCACAATCATACGTCCGCCTCGACGCCTCCGCCCTATCCTTAGAAGGTGCGATGGGGTGCACTAGCGGGGCAAGCGCA

tail rcluster
6758	1232	AATTCGAGTGGTCTCGTAACCGATCGCTCGCTGACTGGGGGTGCATGATGTGCGCAGTACCTAATGCGGCAGTGAGAAAATAAGCCAAAAGCACGCGCATCCTGAAGCGCTAATTG	CGGATGGCAGTAGGCCTGGGCTAGTTGTATCGCATGACTTAATATGATTTAACGGGGACCATCAGCCGGCGGTTTACATATGTAGTTCCTGCGCAAACACTCTTCTGGCTCAGATTTTCACGCCC
6895	1232	AATTCGAGTGGTCACGTAACAGATCGCTCGCTGACTGGGGGTGCATGATGTGCGCAGTACCTAATGCGGCAGTGAGAAAATAAGCCAAGAGCACGCGCATCCTGAAGCGCTAATTG	CGGATGGCAGTAGCCTGGGCTAGTTGTATCGCATGAGTTAATTGATTTAACGGGGACCATTCAGCCGGCGGTTTAACATATGTAGTTCCTGCAGCAACACTCTTTCGGCTCAGATTTTCACGCCC
5906	1233	AATTCATCTGCGTATTCACCCGTGCGAGGCACTTCTACGTAAGGCGTCGACGCAAGCACATAGAAGGTTACGCTGGGTTAAGGCTTGGTTTCTTCTTTAACTTATCCGGGTTATAC	CGGAAGGAGCTACTCGGTTCGGGCGCATCTCTGTATGGTGCCTCGCACTAATCTGTTACCTGTTCTTAGCGAGAGGGCACTCCGGAGTAGCTTAGTCATACAACCGCCCGCTTGCTTCGCCTTTG
6033	1233	AATTCATCTGCGTATTCACCCGTGCGAGGCACTTCTACGTAAGGCGTCGACGCAAGCACATAGAAGGTTACGCTGGGTTAAGGCTTGGTTTCTTCTTTAACTTATCCGGGTTATAC	CGGAAGGAGCTACTCGGTTCGGGCGCATCTCTGTATGGTGCCTCGCACTAATCTGTTACCTGTTATTAGCGAGAGGGCACTCCGGAGTAGCTTAGTCATACAACCGCCCGCTTGCTTCGCCTTTG
7618	1233	AATTCATCTGCGTATTCACCCGTGCGAGGCACTTCTACGTAAGGCGTCGAAGCAAGCACATAGAAGGTTACGCTGGGTTAAGGCTTGGTTTCTTCTTTAACTTATCCGGGTTATAC	CGGAAGGAGCTACTCGGTTCGGGCGCATCTCTGTATGGTGCCTCGCACTAATCTGTTACCTGTTATTAGCGAGAGGGCACTCCGGAGTAGCTTAGTCATACAACCGCCCGCTTGCTTCGCCTTTG
7826	1233	AATTCATCTGCGTATTCACCCGTGCGAGGCACTTCTACGTAAGGCGTCGACGCAAGCACATAGAAGCTTACGCTGGGTTAAGGCTTGGTTTCTTCTTTAACTTATCCGGGTTATAC	CGGAAGGAGCTACTCGGTTCGGGCGCATCTCTGTATGGTGCCTCGCACTAATCTGTTACCTGTTATTAGCGAGAGGGCACTCCGGAGTAGCTTAGTCATACAACCGCCCGCTTGCTTCGCCTTTG
5907	1234	AATTCGCTATCCCAAGGAGCGGTGTGGCGATTCTCGCAGGCTGACCGTGTGACCAAGTCCACGGACGAATTTGCAGGGACGTACGTAGTGACTTATTGGTTTGGTATGGGAGTCGC	CGGTCTGCAGGTTCTATCGTCCCTTTTGCACTCGATTTAAGGGCCTCATTCTTAGATCTTGGCGAAGTCATGCACTGTCGACTGGTGGTACCCGTCACTTATGATCATTGCGCGGAGCGCATCCG
6632	1234	AATTCGCTATCCCAAGGAGCGGTGTGGCGATTCTCGCAGGCTGACCGTGTGACCAAGGCCACGGACGAATTTGCAGGGACGTACGTAGTGACTTATTGGTTTGGTATGGGAGTCGC	CGGTCTGCAGGTTCTATCGTCCCTTTTGCACTCGATTTAAGGGCCTCATTCTTAGATCTTGGCGAAGTCATGCACTGTCGACTGGTGGTACCCGTCACTTATGATCATTGCGCGGAGCGCATCCG
7199	1234	AATTCGCTATCCCAAGGAGCGGTGTGGCAATTCTCGCAGGCGACCGTGTGACCAAGGCCACGGACGAATTCGCAGGGACGTACTAGTGACTATTGGTTTGGTATGGGATCGCCATA	CGGGGTCCAAACTCTGCAGTTCGATCGTCCCTTTTGCACTCGATTTAAGGGCCTCATTTTATGATCTTGGCGAATCTGCCAGCTGTCGACTGGTGGTACCGTCATTCATGGATCTTGCGCGGAGC
7985	1234	AATTCGCTATCCCAAGGAGCGGTGTGGCGATTCTCGCAGGCGACCGTGTGACCAAGGCCACGGACGAATTCGCAGGGACGTACTAGTGACTATTGGTTTGGTATGGGATCGCCATA	CGGGGTCCAAACTCTGCAGTTCGATCGTCCCTTTTGCACTCGATTTAAGGGCCTCATTTTATGATCTTGGCGAATCTGCCAGCTGTCGACTGGTGGTACCGTCATTCATGGATCTTGCGCGGAGC
```

Now we can calculate the number of clusters in rcluster file 

```
cut -f2 rcluster | uniq | wc -l 
1234
```

1234 clusters! 

Now rainbow will split the clusters into smaller clusters representing significant variants. 

```
rainbow div -i rcluster -o rbdiv.out 

less rbdiv.out
1       1       AATTCCAGGGCTGATATACCACCACACAATAGCTGGTTCTGGTCCTCGAATAGTACAATCACCGGAAAGCCAGATACTCCACATTTATTTCATACGCTACTGTGCGACAGCCCTTT    CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC   1       
```

In rainbow, we can use -f flag to control what the min frequency of an allele is necessary to divide it into its own clister. Since this is from multiple individuals, we want to lower this from the default of 0.2.

```
rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
```

This step will matter more when using real data

Now, we can use the PE reads to merge divided clusters to check clustering and dividing of the previous steps, which was based on forward reads.

```
rainbow merge -o rbasm.out -a -i rbdiv.out

head rbasm.out 
E 1
C 0
L 116
S AAAGGGCTGTCGCACAGTAGCGTATGAAATAAATGTGGAGTATCTGGCTTTCCGGTGATTGTACTATTCGAGGACCAGAACCAGCTATTGTGTGGTGGTATATCAGCCCTGGAATT
N 1
R 1:0:0
//
C 1
L 125
S CGGGTTCGTACGGGCGACAATGCTAGGATAGGTTCCTTCCCCGTAAAATGGCCTGCTCGTACCTGAACGTTTGTATGGATATAGGCATCCACTCGGTACGGTCGCCTGCGTCGATCGTTAGAACC
```

We can also use the -r flag, which means min # of reads to assemble. The default is 5, but we are going to use 2 because we working with a reduced dataset!

```
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
```

The rbasm file has both optimal and suboptimal contigs. We can use the following code to extract the optimal contigs. It looks at all the contigs in a cluster. If any of the contigs contain both forward and PE reads, that contig is optimal. If there are no overlapping contigs, the contigs with the most assembled PE reads is output with forward read contig with a 10N spacer. If two contigs have = # of reads, the longer contig is output. 

```
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

If we were using real data, we would stop here and double check rainbow's assembly. But we don't have real data hooray 

rainbow isn't perfect, so now we need to align and cluster by seq similarity the resulting contigs using cd-hit

```
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9

================================================================
Program: CD-HIT, V4.6 (+OpenMP), Jan 11 2018, 15:28:37
Command: cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0
         -T 0 -c 0.9

Started: Tue Mar  1 15:22:27 2022
================================================================
                            Output                              
----------------------------------------------------------------
total number of CPUs in the system is 80
Actual number of CPUs to be used: 80

total seq: 1195
longest and shortest : 371 and 230
Total letters: 300674
Sequences have been sorted

Approximated minimal memory consumption:
Sequence        : 0M
Buffer          : 80 X 12M = 970M
Table           : 2 X 16M = 33M
Miscellaneous   : 4M
Total           : 1009M

Table limit with the given memory limit:
Max number of representatives: 248606
Max number of word counting entries: 10919152

# comparing sequences from          0  to         14
---------- new table with       14 representatives
# comparing sequences from         14  to         28
..............---------- new table with       14 representatives
# comparing sequences from         28  to         42
..............---------- new table with       14 representatives
# comparing sequences from         42  to         56
..............---------- new table with       14 representatives
# comparing sequences from         56  to         69
.............---------- new table with       13 representatives
# comparing sequences from         69  to         82
.............---------- new table with       13 representatives
# comparing sequences from         82  to         95
.............---------- new table with       13 representatives
# comparing sequences from         95  to        108
.............---------- new table with       13 representatives
# comparing sequences from        108  to        121
.............---------- new table with       13 representatives
# comparing sequences from        121  to        134
.............---------- new table with       13 representatives
# comparing sequences from        134  to        146
............---------- new table with       12 representatives
# comparing sequences from        146  to        158
............---------- new table with       12 representatives
# comparing sequences from        158  to        170
............---------- new table with       12 representatives
# comparing sequences from        170  to        182
............---------- new table with       12 representatives
# comparing sequences from        182  to        194
............---------- new table with       12 representatives
# comparing sequences from        194  to        206
............---------- new table with       12 representatives
# comparing sequences from        206  to       1195
....................---------- new table with      989 representatives

     1195  finished       1195  clusters

Apprixmated maximum memory consumption: 1010M
writing new database
writing clustering information
program completed !

Total CPU time 6.31
```

The -c flag sets the % of seq similarlity to group contigs by. In this example, it's 90%. Can also try now with 95%, 85%, 80%, and 99%. 

We can vary the uniq seq copy cutoff and final clustering similarity quickly/easily by using the `remake_reference.sh` script.

Get the script here: 

```
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh
```

Now we can use it to remake the reference and test out new parameters! 

```
bash remake_reference.sh 4 4 0.90 PE 2

remake_reference.sh: line 35: [: an 0: integer expression expected
All required software is installed!
dDocent remake_reference version 2.8.7
dDocent assembled 7989 sequences (after cutoffs) into 1235 contigs
```

Skipping `ReferenceOpt.sh` portion. 

Check out the reference fasta file

```
head reference.fasta
>dDocent_Contig_1
NAATTCATGGGATTCCCTGAGAGCACGAACGTCATTTACATCTAATACTCATTGGCACGTCATGCATTGGCAAAGACGAGTAGTTAGTGATAACGCCTAATCACCACGTCAACTGAANNNNNNNNNNGTGGGGTAGCGGGATAGTAAGCCCCCTGGATTGTCCTATTGACTGCGAAGACAAATAGCAGAGGTTCATACGCTCGGTCCTTTGCGCAGAGGACGGACGATTCGGACGCTATCCTAATTTTGCCGN
>dDocent_Contig_2
NAATTCGTTTGCCCACGGCTTCCACTAAAAGTTGCCCGCAGAACGGATCACTCCAGTATATGCTGCAGTTTGGATGTGAGAGCGGATAGTTTTGCTAATCATCCACGGGCCGTTAGTNNNNNNNNNNCGTAAGGCATCGGATTTACACGGTATACGTCCGATGCATTGCTTGTACCCACGTCCGAATTCATCGACGTGCGCACTCCTGATTATAACTTAACAATAGTCATAACGCCGGGCTCCCTGGTTCCGN
>dDocent_Contig_3
NAATTCACGGCTATCAACTAGGATGGTGGTTACTATTAGTGAGTGCTGTGTATTTCCGCTGCCGTCACTTGCAAGGCAGTAAACCCTTGGTGGCACGTGTAATCCAGCGTATGCAATNNNNNNNNNNCCTGGCCATAGTATTGCTCTAGCATAAAACAAGAGTTATGTATCTTGCCTTCCGGCTAGTCACCTATAGTGATTTGAGCTATTGAAAAGTCACGTTGACTGGAGGTAGAGAGTGGAATACTTCCGN
>dDocent_Contig_4
NAATTCAGCAACTCAAGTAATTCTGTGACTGCCACACCTTTCACCTGTAAGGCACTCGCGTACATCATTAGATCTTATTTGAAAGACCTGGCGTCGCCAATGTTGTCCGCAATAATCNNNNNNNNNNTGCGTACTCACGTTGTTATATAATGCAGCCGTCACACAATTTCGTGGATCGGCTACGGTGCGGGACTGAGACATACGTACGGTCAATAGGAGTAATAATCGCTTCATCATGATACTGGCTGTCCGN
>dDocent_Contig_5
NAATTCTCTGGCATTAATACCTTTATTTCTTTCCCGGAATTTGGTCCATACGCCAAAACCACATTAACTTTACACAGACCATGTCCTCGACGGTCGAATTTAGACAAATTTCTAGTGNNNNNNNNNNCCGTTATTGAACCGAACTATCCTGTTTGATTCGGGGCCTTGGATTTTGACTGGCGTAAGTGCACCGAATTTATAGTATACAATTTTTCACGGGGTAGACGAGTGCGATATCGATCGAGTGAACCGN
```

Do some investigating of the data

```
# number of lines
wc -l reference.fasta
2470 reference.fasta

# number of contigs
zgrep -c ">" reference.fasta
1235
# OR
mawk '/>/' reference.fasta | wc -l 
1235

# check that seqs follow expected format
mawk '/^NAATT.*N*.*CCGN$/' reference.fasta | wc -l
1235
grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l
1235
```

Everything looks good! We now know that there are 1235 loci in the data. 

