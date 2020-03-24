# analysis.sh
# C: Jul  8, 2019
# M: Feb 26, 2020
# A: Leandro Lima <leandro.lima@gladstone.ucsf.edu>


########################################
## Variables with tools and databases ##
########################################


REF_GENOME=/root/resources/chr19.fa
DATA_DIR=/root/data/
GATK_DIR=/root/tools/GenomeAnalysisTK-3.8-1-0


#############
## Mapping ##
#############


cd /root/analysis


# BWA - http://bio-bwa.sourceforge.net/bwa.shtml

bwa mem $REF_GENOME $DATA_DIR/patient3.1.fastq.gz $DATA_DIR/patient3.2.fastq.gz > patient3.sam
bwa mem $REF_GENOME $DATA_DIR/patient4.1.fastq.gz $DATA_DIR/patient4.2.fastq.gz > patient4.sam
bwa mem $REF_GENOME $DATA_DIR/patient5.1.fastq.gz $DATA_DIR/patient5.2.fastq.gz > patient5.sam


# Picard - https://broadinstitute.github.io/picard/command-line-overview.html#Overview

# SAM to BAM, with read group information
for pat_id in patient3 patient4 patient5; do
    java -jar /root/tools/picard.jar AddOrReplaceReadGroups \
          I=$pat_id.sam \
          O=$pat_id.bam \
          RGID=4 \
          RGLB=lib1 \
          RGPL=illumina \
          RGPU=unit1 \
          RGSM=$pat_id 
done


# Samtools - http://www.htslib.org/doc/samtools.html

# Sort
samtools sort patient3.bam -o patient3.sort.bam
samtools sort patient4.bam -o patient4.sort.bam
samtools sort patient5.bam -o patient5.sort.bam

# Create index
samtools index patient3.sort.bam
samtools index patient4.sort.bam
samtools index patient5.sort.bam



# Check file sizes
ls -lh

# Remove SAM file

# Check some flags (number of mapped/unmapped reads, etc.)


#####################
## Variant Calling ##
#####################


# GATK Haplotype Caller - https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller

GATK_DIR=/root/tools/GenomeAnalysisTK-3.8-1-0

for pat_id in patient3 patient4 patient5; do
    java -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R $REF_GENOME \
        -I $pat_id.sort.bam \
        -o $pat_id.gatk.g.vcf \
        --emitRefConfidence GVCF
done


# GATK Joint Genotyping
java -jar $GATK_DIR/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R $REF_GENOME \
    --variant patient3.gatk.g.vcf \
    --variant patient4.gatk.g.vcf \
    --variant patient5.gatk.g.vcf \
    -o all_patients.vcf


################
## Annotation ##
################


# snpEff - http://snpeff.sourceforge.net/SnpEff_manual.html

SNPEFF_DIR=/root/tools/snpEff
# snpEff
java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar \
    -v -stats all_patients.html \
    GRCh38.86 all_patients.vcf > all_patients.ann.vcf

# Adding ID field from dbSNP
java -jar $SNPEFF_DIR/SnpSift.jar annotate \
    -id /root/resources/dbSNP_chr19.vcf.gz \
    all_patients.ann.vcf > all_patients.dbSnp.vcf

rm all_patients.ann.vcf


###################
## APOE analysis ##
###################


# SnpSift extractFields - http://snpeff.sourceforge.net/SnpSift.html#Extract

# Extracting gene name, rsID and genotypes
cat all_patients.dbSnp.vcf | grep CHROM | cut -f1-5,8,10- > APOE_status.txt
java -jar /root/tools/snpEff/SnpSift.jar \
    extractFields \
    all_patients.dbSnp.vcf \
    CHROM POS ID REF ALT EFF[0].GENE GEN[*].GT \
    | grep -E 'rs429358|rs7412' >> APOE_status.txt 

column -t APOE_status.txt
