#!/bin/bash
declare -a StringArray1=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
for fq in "${StringArray1[@]}"
do
gzip -c $fq*.assembled.fastq > $fq.norm.merge.assembled.fq.gz
done