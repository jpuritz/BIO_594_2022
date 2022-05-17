#!/bin/bash
for i in "${StringArray[@]}"
do
fastp --detect_adapter_for_pe --thread 10 -V -e 30 -q 30 -p --html /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.html --json /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}.json -i /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-1_2-1/${i}.F.fq.gz -I /home/jgreen/repos/BIO594_work/course_project/data/gDNA/demux/capture_1-1_2-1/${i}.R.fq.gz -o /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.F.fq.gz -O /home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/${i}_trim.R.fq.gz
done