#!/bin/bash
#BUSCO
busco -i trinity.transcripts.cluster90.fasta -o mollusc_oases_busco -l ~/databases/mollusca_odb10/ -m transcriptome --cpu 5

#Detonate
rsem-eval-estimate-transcript-length-distribution trinity_out_dir.Trinity.fasta trinity.rsem.eval.ETLD.param
rsem-eval-calculate-score -p 5 --transcript-length-parameters trinity.rsem.eval.ETLD.param --paired-end cDNA.F.fq.gz cDNA.R.fq.gz trinity.transcripts.cluster90.fasta trinity_cluster 638

#Transrate
transrate --assembly trinity_out_dir.Trinity.fasta --left cDNA.F.fq.gz --right cDNA.R.fq.gz  --threads 5 

#RNAquast
rnaQUAST.py -c trinity_out_dir.Trinity.fasta -1 cDNA.F.fq.gz -2 cDNA.R.fq.gz --busco ~/databases/mollusca_odb10/ -t 10 -o trinity_rnaquast
