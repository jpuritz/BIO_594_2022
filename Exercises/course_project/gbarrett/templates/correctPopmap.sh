#!/bin/bash

bcftools query -l $vcf > vcfsamples.txt

awk -F',' 'NR==FNR{a[\$1];next} \$4 in a {print \$5 \$9}' vcfsamples.txt $meta > popmap.tsv

# revise populations For BayPass group

#sed -E -e 's/IMR1|IMO/Itk_high/' -e 's/IMR2|BJANE|LIMS/Itk_highM/' -e 's/OksR1|IMR3|IM3|LL3/Itk_lowM/' -e 's/Itk4|Itk3|Itk4.5/Itk_low/' \
#-Ee 's/GCL|CGK|Kup2/Kup_high/' -e 's/LTER|KupR2//' popmap.tsv