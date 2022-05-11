#!/bin/bash
# Subset Individuals
#########################
vcftools --vcf "${vcf}" --remove "${probSamples}" \
--recode --recode-INFO-all --stdout > fb.F.vcf

bcftools query -l "fb.F.vcf" > "final_indv.txt"

# based on watershed
awk -F',' '{print \$1"\t"\$3}' "$meta" > "watershed.txt"

# Site-Based Filters
########################
vcftools --vcf fb.F.vcf --mac 6 --minQ 20 \
--minDP 3 --min-meanDP 3 \
--maxDP 67 --max-meanDP 67 \
--recode --recode-INFO-all --stdout > fb.F.mac.6.minQ20.minDP3.maxDP67.vcf

# PER-POP MISSINGNESS
########################## 
sh ${baseDir}/bin/filter_miss_by_pop.sh fb.F.mac.6.minQ20.minDP3.maxDP67.vcf "watershed.txt" .65 1 > fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.vcf

# INFO
########
bcftools filter -i'AB>0.25 && AB<0.75 | AB < 0.01' fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.vcf | \
bcftools filter -i'SAF / SAR > 100 && SRF / SRR > 100 | SAR / SAF > 100 && SRR / SRF > 100' | \
bcftools filter -i'MQM / MQMR > 0.9 && MQM / MQMR < 1.05' | \
bcftools filter -i'PAIRED > 0.05 && PAIREDR > 0.05 && PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05' | \
bcftools filter -i'QUAL/ INFO/DP > .25' > fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.vcf

# HWE !Decompose into SNPs!
############################
vcfallelicprimitives --keep-info --keep-geno fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.vcf > fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.snps.vcf
vcftools --vcf fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.snps.vcf --remove-indels --recode --recode-INFO-all --stdout > fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.snps.ri.vcf

# per-pop hardy-weinberg eq.
##############################
vcf_final=fb.F.mac.6.minQ20.minDP3.maxDP67.miss65.info.snps.ri.vcf
rm badloci
pops="\$(awk -F , 'FNR>1 {print \$1}' $meta | sort -u)"
# adding double quotes prints every element in one line
for id in \${pops[@]}; do
    awk -F, -v pop="\$id" '{ if (\$1 == pop) print \$5}' ${meta} > "\${id}.indv"
    vcftools --vcf \${vcf_final} --keep "\$id.indv" --hardy --out \$id
    awk -v x=.0001 '{if(\$6 < x) print \$1 \$2} ' \$id.hwe >> badloci
    loc="\$(cat badloci | wc -l)"
done
cat badloci | sort -u > loci.to.remove
vcftools --vcf \${vcf_final} --exclude-positions loci.to.remove --recode --recode-INFO-all --stdout > consensus.vcf




