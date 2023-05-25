#!/bin/bash

# This script should be run on your local laptop or computer
# Please make sure you have installed Docker Desktop using the instructions at: https://www.docker.com/products/docker-desktop/
# Before running this script, do the following
##  1. open the Docker Dekstop and make sure it is running
##	   If Docker Desktop gives a pop-up that it requires a newer WSL kernel version, open a new command prompt window and run the below command:
##	   $ wsl --update
##	   Restart Docker Desktop
##  2. Download the workshop materials at https://github.com/gladstone-institutes/Bioinformatics-Workshops/raw/master/intro-rna-seq/Intro_to_RNA-seq_data_analysis.zip?raw=true.
##  3. Unzip the workshop materials in the Downloads folder
##  4. Open a new command prompt window 

##run the below commands in the command prompt window
#check if docker is running
#the below command will print out the docker version if docker is running
docker --version

#get the docker image from https://hub.docker.com/r/nfcore/rnaseq/
#we will be using this image for all the analyses
docker pull nfcore/rnaseq 

#check if the docker image was downloaded successfully
docker images

#spin up a docker container
#change the path below to the path on your computer where the downloaded materials are unzipped
docker run --name rna_bash --rm -it -v C:\Users\ayushi.agrawal\Downloads\Intro_to_RNA-seq_data_analysis\Intro_to_RNA-seq_data_analysis:/home nfcore/rnaseq bash

#a new prompt will apear and the conatiner is now active
#the home directory should have all the contents of the downloaded materials folder
#go to the home directory in the docker container
cd home

#run fastqc
fastqc Bacteria_GATTACA_L001_R1_001.fastq

#trim the reads using cutadapt
cutadapt -a file:Adapter_Sequence.fasta -o trimmed.fastq Bacteria_GATTACA_L001_R1_001.fastq

#run fastqc on the trimmed reads
fastqc trimmed.fastq


#STAR does not work on windows so, we will be using another aligner "HISAT2" for the demo
#create a new folder
mkdir hisat2_index

#create the hisat2 index
hisat2-build -p 7 rDNA_sequence.fasta hisat2_index/rDNA_

#run hisat2 for the trimmed reads
#HISAT2 and STAR are different aligners with different defaults, scoring algorithms, etc. This might result in different outputs from the two aligners. 
hisat2 -p 7 -x hisat2_index/rDNA_ -U ./trimmed.fastq -S hisat2_aligned_out.sam

#generate the read count matrix using featureCounts
featureCounts -a rDNA.gtf -t CDS -o counts.txt hisat2_aligned_out.sam


# END #