#!/bin/bash

# This script should be run on your local laptop or computer
# Please make sure you have installed Docker Desktop using the instructions at: https://www.docker.com/products/docker-desktop/
# Before running this script, do the following
##  1. open the Docker Dekstop and make sure it is running
##  2. open a new terminal (MacOS) or command prompt (Windows) window 

#check if docker is running
#the below command will print out the docker version if docker is running
docker --version

#get the docker image from https://hub.docker.com/r/nfcore/rnaseq/
#we will be using this image for all the analyses
docker pull nfcore/rnaseq 

#check if the docker image was downloaded successfully
docker images

#spin up a docker container
##for M1 chip mac, use the below command
docker run --platform linux/amd64 --name rna_bash --rm -it -v ~/Downloads/Intro_to_RNA-seq_data_analysis:/home nfcore/rnaseq bash
##for all other computers, use the below command
docker run --name rna_bash --rm -it -v ~/Downloads/Intro_to_RNA-seq_data_analysis:/home nfcore/rnaseq bash

#a new prompt will apear and the conatiner is now active
#the home directory should have all the contents of ~/Downloads/Intro_to_RNA-seq_data_analysis
#go to the home directory in the docker container
cd home

#run fastqc
fastqc Bacteria_GATTACA_L001_R1_001.fastq

#trim the reads using cutadapt
cutadapt -a file:Adapter_Sequence.fasta -o trimmed.fastq Bacteria_GATTACA_L001_R1_001.fastq

#run fastqc on the trimmed reads
fastqc trimmed.fastq

#create a new folder
mkdir star_index

#create the STAR index
STAR --runMode genomeGenerate --genomeDir ./star_index --genomeFastaFiles rDNA_sequence.fasta --genomeSAindexNbases 3

#run STAR for the trimmed reads
STAR --genomeDir ./star_index --readFilesIn ./trimmed.fastq

#generate the read count matrix using featureCounts
featureCounts -a rDNA.gtf -t CDS -o counts.txt Aligned.out.sam




