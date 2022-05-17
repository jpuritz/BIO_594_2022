#!/bin/bash
declare -a StringArray=("CASE" "CA" "CON" "SE")
cd-hit-est -i "$i".cDNA.assembled.fastq -o "$i".cDNA.c95.cluster -c 0.95 -n 10 -d 0 -M 16000 -T 8
