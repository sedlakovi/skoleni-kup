#!/bin/bash

# Specify sample name
SAMPLE_NAME=test

# Set output directories for various data types
FASTQ_DIR=~/str_project/data/reads
TRIMMED_FASTQ_DIR=~/str_project/data/trimmed_reads
BAM_DIR=~/str_project/data/bam
VCF_DIR=~/str_project/data/vcf
STR_DIR=~/str_project/data/str_calls

# Create those directories if they do not exist
cd
mkdir -p $FASTQ_DIR
mkdir -p $TRIMMED_FASTQ_DIR
mkdir -p $BAM_DIR
mkdir -p $VCF_DIR


# Set path to reference files
REFERENCE=/srv/resources/hg19.fa.gz
CODIS_BED=/srv/resources/ref_test.bed
DBSNP_BED=/srv/resources/dbSnp155Common.bed

# Step 1. Trim FASTQ files
fastp -i "$FASTQ_DIR"/"$SAMPLE_NAME"_R1.fastq -I "$FASTQ_DIR"/"$SAMPLE_NAME"_R2.fastq -o "$TRIMMED_FASTQ_DIR"/"$SAMPLE_NAME"_trimmed_R1.fastq -O "$TRIMMED_FASTQ_DIR"/"$SAMPLE_NAME"_trimmed_R2.fastq

# Step 2. Align trimmed reads to the reference genome
bwa mem -t 4 $REFERENCE "$TRIMMED_FASTQ_DIR"/"$SAMPLE_NAME"_trimmed_R1.fastq  \
"$TRIMMED_FASTQ_DIR"/"$SAMPLE_NAME"_trimmed_R2.fastq | samtools view -b | samtools sort > "$BAM_DIR"/"$SAMPLE_NAME".bam

# Step 3. Variant calling - SNPs and INDELs
bcftools mpileup -a AD,DP,SP -Ou -f $REFERENCE \
"$BAM_DIR"/"$SAMPLE_NAME".bam | bcftools call -f GQ,GP \
-mO z -o "$VCF_DIR"/"$SAMPLE_NAME".vcf.gz

# Step 4. Annotate VCF files - add rsID
bcftools annotate -c CHROM,FROM,TO,ID -a $DBSNP_BED  -o "$VCF_DIR"/"$SAMPLE_NAME"_annotated.vcf.gz "$VCF_DIR"/"$SAMPLE_NAME".vcf.gz

# Step 5. Variant calling - STR
samtools addreplacerg -r "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tLB:lib1" "$BAM_DIR"/"$SAMPLE_NAME".bam -o "$BAM_DIR"/"$SAMPLE_NAME"_rg.bam
samtools index "$BAM_DIR"/"$SAMPLE_NAME"_rg.bam
HipSTR --bams "$BAM_DIR"/"$SAMPLE_NAME"_rg.bam --fasta $REFERENCE --regions $CODIS_BED --str-vcf "$STR_DIR"/"$SAMPLE_NAME".vcf.gz --max-str-len 200 --min-reads 1