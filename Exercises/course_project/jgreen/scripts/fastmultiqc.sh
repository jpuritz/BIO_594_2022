#!/bin/bash
fastqc -f fastq *.fq.gz
multiqc -o multiqc .