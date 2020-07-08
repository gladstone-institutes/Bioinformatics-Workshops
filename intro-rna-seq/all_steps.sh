#!/bin/bash

fastqc Bacteria_GATTACA_L001_R1_001.fastq

cutadapt \
-a file:Adapter_Sequence.fasta \
-o trimmed.fastq \
Bacteria_GATTACA_L001_R1_001.fastq

fastqc trimmed.fastq

mkdir star_index

STAR \
--runMode genomeGenerate \
--genomeDir ./star_index \
--genomeFastaFiles rDNA_sequence.fasta \
--genomeSAindexNbases 3

STAR \
--genomeDir ./star_index \
--readFilesIn ./trimmed.fastq

featureCounts \
-a rDNA.gtf \
-t CDS \
-o counts.txt \
Aligned.out.sam
