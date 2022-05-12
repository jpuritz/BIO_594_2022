#!/bin/bash
for fq in *assembled.fastq
do
echo $fq
cat $fq | echo $((`wc -l`/4))
done