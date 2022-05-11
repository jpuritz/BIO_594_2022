#!/bin/bash
env="\$(basename ${lfmm} .vcf | sed -e 's,lfmm_,,' -e 's,.env,,')"

# Remove NA's
awk '!/NA/ {print \$1}' ${lfmm} > "\${env}.indv"

# write env_vcf file
bcftools view -I --force-samples -S "\${env}.indv" ${vcf} > "\${env}.vcf" # force samples not present in vcf header, -I = don't recaculate INFO fields

# Get FINAL individual list 
bcftools query -l "\${env}.vcf" > "\${env}_${vcf.baseName}.indv"

# write ped file
plink --vcf "\${env}.vcf" --recode --allow-extra-chr --out "\${env}"

# Remove absent individuals in group_list from FINAL
awk 'NR==FNR{a[\$1];next} \$1 in a' "\${env}_${vcf.baseName}.indv" $group > "\${env}.grp"
# convert vcf to baypass script [group] [output prefix]
cat "\${env}.vcf" | perl ${baseDir}/bin/vcf2baypass.pl "\${env}.grp" "\${env}.frq"

# Remove absent individuals in LEA.env from FINAL
awk -F'[ ]' 'NR==FNR{a[\$1];next} \$1 in a {print \$2}' "\${env}_${vcf.baseName}.indv" $lfmm > "\${env}.env"



