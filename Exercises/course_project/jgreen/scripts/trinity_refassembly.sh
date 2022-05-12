#!/bin/bash
## Mapping gDNA reads to de novo transcriptome

bwa index -a bwtsw trinity.transcripts.cluster90.fasta trinity.bwa
bwa mem -t 20 trinity.transcripts.cluster90.fasta gDNA.F.fq gDNA.R.fq | samtools view -bu - | samtools sort -@4 - -o trinity.c90.sorted.bam

## Trinity transcriptome guided assembly

Trinity --genome_guided_bam trinity.c90.sorted.bam \
         --genome_guided_max_intron 5000 \
         --max_memory 50G --CPU 20


~/miniconda3/envs/trinity/util/misc/process_GMAP_alignments_gff3_chimeras_ok.pl \
     --genome trinity_out_dir.Trinity.fasta \
     --transcripts trinity_out_dir/<filename>.fasta \
     --SAM | samtools view -Sb | samtools sort -o trinity-GG.gmap.bam
