#!/bin/bash

#This script should be run on the dev3 node of the UCSF Wynton HPC cluster
#before running this script, do the following
##  1. transfer the workshop files from your local computer to your wynton account
##     {local}$ scp -r Downloads/Intro_to_RNA-seq_data_analysis/ alice@dt2.wynton.ucsf.edu:~
##  2. create the singularity container rna_seq_container.sif on wynton using the below command
##     [alice@log1 ~]$ singularity build rna_seq_container.sif docker://nfcore/rnaseq 

cd ~/Intro_to_RNA-seq_data_analysis/

singularity exec rna_seq_container.sif fastqc Bacteria_GATTACA_L001_R1_001.fastq

singularity exec rna_seq_container.sif cutadapt \
-a file:Adapter_Sequence.fasta \
-o trimmed.fastq \
Bacteria_GATTACA_L001_R1_001.fastq

singularity exec rna_seq_container.sif fastqc trimmed.fastq

singularity exec rna_seq_container.sif mkdir star_index

singularity exec rna_seq_container.sif STAR \
--runMode genomeGenerate \
--genomeDir ./star_index \
--genomeFastaFiles rDNA_sequence.fasta \
--genomeSAindexNbases 3

singularity exec rna_seq_container.sif STAR \
--genomeDir ./star_index \
--readFilesIn ./trimmed.fastq

singularity exec rna_seq_container.sif featureCounts \
-a rDNA.gtf \
-t CDS \
-o counts.txt \
Aligned.out.sam
