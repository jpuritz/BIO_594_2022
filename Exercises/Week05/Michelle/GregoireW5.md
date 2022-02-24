# How to determine how many loci are in a RadSeq dataset using dDocent
Loci are the genetic location of genes on chromosomes. In population genetics, loci can be used to see differences among populations. 

RADSeq (Restriction site Associated DNA Sequencing) is a technique that is often used to sequence DNA of populations in a relatively cheap manner. 
Analysis of RADSeq data can be performed using the dDocent pipeline. Below you will find how to use dDocent commands to determine the number of contigs that designate the loci assembled by dDocent. These commands are to be run on files that have already been demultiplexed and named appropirately.

## Code to download dDocent bash scripts
> curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh
> curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ReferenceOpt.sh
The first line of code will download the script required to build the reference genome and the second will download a script that will assemble references from an interval of cutoff values and percent sequence similarity values.

## Codes to run dDocent for reference assembly
> bash remake_reference.sh 10 10 0.95 PE 2
Here the code uses cutoff values of 10 and the paired-end assembly algorithm (PE --this is for RADseq data that does not have random shearing ie: ddRAD, ezRAD)
This will output a file called reference.fasta. The number of assembled contigs (or loci) can be determined by using an AWK command and a wc -l flag as follows
> mawk '/>/' reference.fasta | wc - l
The following code can be used to optimize the parameters for the assembly if desired:
> bash ReferenceOpt.sh 4 8 4 8 PE 16
The numbers after the script name and before specifying the paired ends indicate the range of cutoff values for the reads.

# How many loci were detected in Week5 take-home exercise data? 
With the parameters that I ran in the code above, 1234 contigs were assembled. These are the unique loci that were sequenced and assembled. 


