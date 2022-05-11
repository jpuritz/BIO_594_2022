#!/bin/bash


awk -F, '{for(i=9;i<=NF-1;i++) printf $i" "; print ""}' $1



#awk -F, '
#NR==1 {
#    for (i=1; i<=NF; i++) {
#        f[$i] = i
#    }
#}
#{ print $(f["bioinformatics_id"]), 
#$(f["H3_JulyAug_Temp_Max"]),
#$(f["max_dist"]),
#$(f["location"])
#}' $1 #| sed 's/  / NA /g' 


for i in *.grp; do
    
    cat "\${env}.vcf" | perl ${baseDir}/bin/vcf2baypass.pl "\${env}.grp" "\${env}.frq"



done
