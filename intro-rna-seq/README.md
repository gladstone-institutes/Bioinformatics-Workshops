# intro-rna-seq
[Link to wiki page](https://github.com/gladstone-institutes/Bioinformatics-Workshops/wiki/Introduction-to-RNA-Seq-Analysis)

### Description of files
1. Single_read.fastq (fastq file with a single read to understand the fastq file format)
2. Bacteria_GATTACA_L001_R1_001.fastq (single-end small practice data with 100k reads for the demo in the workshop)
3. Adapter_Sequence.fasta (fasta file with adapter sequence for demo with cutadapt)
4. rDNA_sequence.fasta (fasta file with the reference genome sequence for demo with STAR aligner)
5. rDNA.gtf (GTF file with the annotations for demo with featureCounts)
6. all_steps_wynton.sh (shell script for running all the analysis steps on UCSF Wynton  command-line interface using the practice data provided)
7. steps_on_wynton_part1.txt (text file with steps used on wynton to setup the folders, upload the data and create a singularity container)
8. steps_on_wynton_part2.txt (text file with steps used on wynton to run the bulk RNA-seq analysis using the demo files)
9. all_steps_docker_desktop.sh (shell script with commands for running all the analysis steps using Docker Desktop and the demo files)
