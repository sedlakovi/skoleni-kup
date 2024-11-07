#!/bin/bash

set -exuo pipefail

apt-get update
apt-get install -y \
    libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2t64 libxi6 libxtst6 \
    makepasswd unzip \
    bedtools bwa tree samtools bcftools fastqc \
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

# Install Trimmomatic
cd /srv
curl -LO https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip
echo 'alias trimmomatic="java -jar /srv/Trimmomatic-0.39/trimmomatic-0.39.jar"' > /etc/profile.d/trimmomatic.sh

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

# Download example FASTQ
mkdir fastq
cd fastq
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz

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
