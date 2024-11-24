#!/bin/bash

apt-get install -y \
    libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2t64 libxi6 libxtst6 \
    unzip \
    bedtools bwa tree samtools bcftools tabix fastp \
    byobu

# Download regex example data
curl -L https://oc.embl.de/index.php/s/t3V7FR21N8uuV9L/download -o regex_data.zip
unzip regex_data.zip
rm regex_data.zip
mv example_files/ regex_data

# Download CLI example data
curl -LO https://swcarpentry.github.io/shell-novice/data/shell-lesson-data.zip
unzip shell-lesson-data.zip
rm shell-lesson-data.zip
mv shell-lesson-data data-shell

# Download reference hg19 FASTA
sudo mkdir -p /srv/resources
cd /srv/resources
sudo curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
sudo gunzip hg19.fa.gz
sudo bgzip hg19.fa.gz

# Create BWA index
sudo bwa index /srv/resources/hg19.fa.gz
sudo samtools faidx /srv/resources/hg19.fa.gz

# Download bigBed to annotate VCF files with rsID
sudo curl -O https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
sudo chmod +x bigBedToBed
sudo curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155Common.bb
sudo ./bigBedToBed dbSnp155Common.bb dbSnp155Common.bed
sudo rm dbSnp155Common.bb

# Download CODIS bed file
sudo wget https://github.com/AnJingwd/STRsearch/raw/refs/heads/master/example/ref_test.bed
sudo sed -i '1d' ref_test.bed

# Download example FASTQ
cd ~
mkdir str_project
cd str_project
mkdir data
cd data
mkdir reads
cd reads
curl -O https://raw.githubusercontent.com/AnJingwd/STRsearch/refs/heads/master/example/test_data/test_R2.fastq
curl -O https://raw.githubusercontent.com/AnJingwd/STRsearch/refs/heads/master/example/test_data/test_R1.fastq

cd ~