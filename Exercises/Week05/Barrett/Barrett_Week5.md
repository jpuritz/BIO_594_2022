# Take-home Exercise for Week 5
My method to obtain the number of loci present in the suite of files is as follows:

```
#!/usr/bin/env bash
PROG=/usr/local/bin
DATA=/home/BIO594/DATA/Week5
RBDIR=/home/gbarrett/BIO_594_2022/Exercises/Week05/Barrett
mkdir -p $RBDIR
for ind in $DATA/*F.fq.gz;
do
	FQ1=${ind}
	FQ2=$(echo $FQ1 | sed -e 's/F.fq.gz/R.fq.gz/')
    echo "Read 1: $FQ1 Read 2: $FQ2"
    rainbow cluster -1 $FQ1 -2 $FQ2 -m 4 -e 2000 > $RBDIR/rbcluster.out 2> $RBDIR/rbcluster.log
done
```
**Cluster** reads with a given number of mismatches into a one group
```
rainbow div -i $RBDIR/rbcluster.out -o $RBDIR/rbdiv.out -k 7 -K 10 -f .5 2> $RBDIR/rbdiv.log
```
**Recursively divide** calculate the count of minor base K and its frequency F and detect minor bases via max(K) / F are lower than the minimum requirement. 
```
rainbow merge -i $RBDIR/rbdiv.out -o $RBDIR/rbmerge.out -a -r 2 -N1000 -R1000 -l 20 -f .9 2> $RBDIR/rbmerge.log
```
**Merge** into clusters: two reads or contigs with the largest overlap are chosen and merged into one consensus sequence
```
# Reference, select best contigs, index
REF=$RBDIR/rainbow.fasta
$PROG/select_best_rbcontig_plus_read1.pl $RBDIR/rbmerge.out $RBDIR/rbdiv.out > $REF
$PROG/samtools faidx ${RBDIR}/rainbow.fasta -o ${RBDIR}/rainbow
INDEX=$RBDIR/rainbow
bwa index -a bwtsw $REF -p $INDEX
```
**Select** the longest contigs as the representative reference and index the reference for downstream analysis

## How many loci are there?

Using `grep -c '<' rainbow.fasta` I obtained **1230** loci. Although, it would appear that 1235 is the correct number. I attribute discrepancies to parameters specified through the workflow. 