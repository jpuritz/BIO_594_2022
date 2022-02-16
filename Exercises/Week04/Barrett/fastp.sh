#!/usr/bin/env bash
INDIR=/home/BIO594/DATA/Week4
OUTDIR=/home/gbarrett/BIO_594_2022/Exercises/Week04/Barrett
#mkdir -p Barrett
for i in $INDIR/*F.fq.gz;
do
	# Input Reads
	fq1=${i}
	fq2=$(echo ${i} | sed 's/F./R./')
    echo $fq1 $fq2
	# Output Reads
	FQ1=$(echo ${fq1} | sed -e 's/F./out_F./' -e "s+$INDIR+$OUTDIR+")
	FQ2=$(echo ${fq2} | sed -e 's/R./out_R./' -e "s+$INDIR+$OUTDIR+")	
    echo $FQ1 $FQ2
	# Output Reports
	OUTPUT=$(basename $fq1 .F.fq.gz)
	echo $OUTPUT
    fastp -i $fq1 -I $fq2 -o $FQ1 -O $FQ2 \
    --cut_right --cut_mean_quality=25 --cut_window_size 6 \
	--detect_adapter_for_pe \
	--correction -h ${OUTDIR}/$OUTPUT.html &> ${OUTDIR}/$OUTPUT.log
done