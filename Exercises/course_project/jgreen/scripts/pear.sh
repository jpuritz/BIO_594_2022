#!/bin/bash
base="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/merge/"
declare -a StringArray=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
for i in "${StringArray[@]}"
do
pear -j 20 -f "$i"_trim.F.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq -r "$i"_trim.R.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq -o "$base""$i"_trim_norm_merge.fq.gz
done