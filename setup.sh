#!/bin/bash

set -exuo pipefail

apt-get update
apt-get install -y \
    libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2t64 libxi6 libxtst6 \
    makepasswd unzip \
    bedtools bwa tree samtools bcftools fastqc tabix fastp \
    ttyd byobu ssh-import-id certbot nginx python3-certbot-nginx

byobu-enable
# manually: byobu

adduser krab --gecos "" --disabled-password
adduser krab sudo --gecos "" --disabled-password

adduser anastazie
adduser anastazie sudo

sed -i 's/^%sudo.*/%sudo   ALL=\(ALL:ALL\) NOPASSWD:ALL/' /etc/sudoers

sudo -iu anastazie bash -c 'ssh-import-id gh:anastazie'
sudo -iu krab bash -c 'ssh-import-id gh:crabhi'


# Download snakemake example data
git clone https://github.com/snakemake/snakemake-tutorial-data.git

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
sudo mkdir data
cd data
sudo mkdir reads
cd reads
sudo curl -O https://raw.githubusercontent.com/AnJingwd/STRsearch/refs/heads/master/example/test_data/test_R2.fastq
sudo curl -O https://raw.githubusercontent.com/AnJingwd/STRsearch/refs/heads/master/example/test_data/test_R1.fastq

cd ~

# Install conda
CONDA_VERSION=2024.10-1
curl https://repo.anaconda.com/archive/Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh -o /srv/install-anaconda.sh

USERS=(adeladanielova ivakulichova janmasek2 janacechova2 jananovackova luboskanca luciecuchalova martinasekowska martinavalisova michalberan2 pavelcapek vlastimilstenzl zbynekdolejsi)

set +x
for u in "${USERS[@]}"; do
    PASS=$(makepasswd --chars=20)
    adduser $u --gecos "" --disabled-password > /dev/null 2>&1
    echo -e "${PASS}\n${PASS}\n" | passwd $u > /dev/null 2>&1

    echo -e "$u\t$PASS"
done

echo
echo "-------------------------"
echo

set -x

for u in krab anastazie "${USERS[@]}"; do
    sudo -u "$u" bash /srv/install-anaconda.sh -b -p /home/$u/anaconda3
    su -l "$u" -c bash -c 'eval "$($HOME/anaconda3/bin/conda shell.bash hook)" && conda init'
done


certbot --agree-tos  -d skoleni-kup.sedlakovi.org --nginx --no-eff-mail -m filip+skolenikup@sedlakovi.org

cp nginx-site-default /etc/nginx/sites-available/default
systemctl restart nginx

cat > /etc/default/ttyd <<EOF
# /etc/default/ttyd

TTYD_OPTIONS="-i lo -W -p 7681 -O login"
EOF
systemctl restart ttyd
