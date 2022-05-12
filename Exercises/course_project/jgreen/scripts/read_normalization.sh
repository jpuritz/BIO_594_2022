#!/bin/bash
#base="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/norm/"
#reads="/home/jgreen/repos/BIO594_work/course_project/data/cDNA/trim/"
#baseg="/home/jgreen/repos/BIO594_work/course_project/data/gDNA/norm/"
#readsg="/home/jgreen/repos/BIO594_work/course_project/data/gDNA/trim/"
#declare -a StringArray1=("CA_J06" "CA_J08" "CA_J11" "CA_J18" "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CON_J02" "CON_J05" "CON_J10" "SE_J01" "SE_J04" "SE_J07")
#declare -a StringArray=("CA_J06_trim" "CA_J08_trim" "CA_J11_trim" "CA_J18G_trim" "CA_J18_trim" "CASE_J03_trim" "CASE_J09_trim" "CASE_J12_trim" "CASE_J13G_trim" "CASE_J13_trim" "CON_J02_trim" "CON_J05_trim" "CON_J10_trim" "CON_J15_trim" "IS_01_trim" "IS_02_trim" "IS_03_trim" "IS_04_trim" "SE_J01_trim" "SE_J04_trim" "SE_J07_trim" "SE_J14G_trim" "SE_J14_trim")
for i in "${StringArray[@]}"
do
insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 30 --left "$i".F.fq.gz --right "$i".R.fq.gz --pairs_together --PARALLEL_STATS --CPU 10 --output "$baseg"
done