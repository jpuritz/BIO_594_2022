#!/bin/bash
for ((i=1;i<=8;i++));
do
process_shortreads -1 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool$i/cDNApool-${i}_R1_001.fastq.gz -2 /home/jgreen/repos/BIO594_work/course_project/data/cDNA/pool$i/cDNApool-${i}_R2_001.fastq.gz -b /home/jgreen/repos/BIO594_work/course_project/data/cDNA/barcodes.txt --inline_inline -D --barcode_dist_1 2 --barcode_dist_2 2 -o /home/jgreen/repos/BIO594_work/course_project/data/cDNA/demux/pool$i/
done