#!/bin/sh -eux

mkdir -p prerequisites/ bin/ java/

#############
cd prerequisites/

wget https://netactuate.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar -xjvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa ../../bin/
cd -

wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar -xjvf samtools-1.11.tar.bz2
cd samtools-1.11
make
cp samtools ../../bin/
cd -

wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
tar -xjvf  bcftools-1.11.tar.bz2
cd bcftools-1.11
make
cp bcftools ../../bin/
cd -

wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
tar -xjvf htslib-1.11.tar.bz2
cd htslib-1.11
make
cp tabix ../../bin/
cd -

wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
tar -xzvf samblaster-v.0.1.26.tar.gz
cd samblaster-v.0.1.26
make
cp samblaster ../../bin/
cd -

wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
tar -xzvf vcftools-0.1.16.tar.gz
cd vcftools-0.1.16
./configure
make
cp src/cpp/vcftools ../../bin/
cd -

wget https://github.com/nferraz/st/archive/v1.1.4.tar.gz -O st-1.1.4.tar.gz
tar -xzvf st-1.1.4.tar.gz
cd st-1.1.4
perl ./Makefile.PL INSTALL_BASE=../../ ; make ; make test ; make install
cd -

################### 

wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
cp picard.jar ../java/

wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
unzip gatk-4.1.9.0.zip
cp gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar ../java/gatk.jar

wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip
unzip haplogrep.zip
cp haplogrep.jar ../java/

wget https://github.com/seppinho/mutserve/releases/download/v1.3.4/mutserve-1.3.4.jar
cp mutserve-1.3.4.jar ../java/mutserve.jar
cd ../

#################

cd RefSeq/
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O hs38DH.fa
../bin/samtools faidx hs38DH.fa
cd -
