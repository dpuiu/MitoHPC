#!/bin/sh -eux

SDIR=`dirname $0`	 # script directory
source $SDIR/init.sh
mkdir -p $HDIR/prerequisites/ $BDIR/ $JDIR/

#############
cd $HDIR/prerequisites/

if [ ! -s $BDIR/bwa ] ; then
  wget https://netactuate.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
  tar -xjvf bwa-0.7.17.tar.bz2
  cd bwa-0.7.17
  make ; cp bwa $BDIR/
  cd -
fi

if [ ! -s $BDIR/samtools ] ; then
  wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
  tar -xjvf samtools-1.11.tar.bz2
  cd samtools-1.11
  ./configure --prefix=$HDIR/ ; make ; make test; make install
  cd -
fi

if [ ! -s $BDIR/bcftools ] ; then
  wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
  tar -xjvf  bcftools-1.11.tar.bz2
  cd bcftools-1.11
  ./configure --prefix=$HDIR/ ; make ; make test; make install
  cd -
fi

if [ ! -s $BDIR/tabix ] ; then
  wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
  tar -xjvf htslib-1.11.tar.bz2
  cd htslib-1.11
  ./configure --prefix=$HDIR/ ; make ; make test; make install
  cd -
fi

if [ ! -s $BDIR/samblaster ] ; then
  wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
  tar -xzvf samblaster-v.0.1.26.tar.gz
  cd samblaster-v.0.1.26
  make ; cp samblaster $BDIR/
  cd -
fi

if [ ! -s $BDIR/vcftools ] ; then
  wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
  tar -xzvf vcftools-0.1.16.tar.gz
  cd vcftools-0.1.16
  ./configure --prefix=$HDIR/ ; make ; make install
  cd -
fi

if [ ! -s $BDIR/st ] ; then
  wget https://github.com/nferraz/st/archive/v1.1.4.tar.gz -O st-1.1.4.tar.gz
  tar -xzvf st-1.1.4.tar.gz
  cd st-1.1.4
  perl ./Makefile.PL INSTALL_BASE=$HDIR/ ; make ; make test ; make install
  cd -
fi

if [ ! -s $BDIR/bedtools ] ; then
  wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
  tar -xzvf bedtools-2.29.2.tar.gz
  cd bedtools2/
  make install prefix=$HDIR/
  cd -
fi

if [ ! -s $BDIR/fastp ] ; then
  wget https://github.com/OpenGene/fastp/archive/v0.20.1.tar.gz -O fastp-0.20.1.tar.gz
  tar -xzvf fastp-0.20.1.tar.gz
  cd fastp-0.20.1/
  make ; cp fastp $BDIR/
  cd -
fi

################### 
if [ ! -s $JDIR/picard.jar ] ; then
  wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
  cp picard.jar $JDIR/
fi

if [ ! -s $JDIR/gatk.jar ] ; then
  wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
  unzip gatk-4.1.9.0.zip
  cp gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar $JDIR/gatk.jar
fi

if [ ! -s $JDIR/haplogrep.jar ] ; then
  wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip
  unzip haplogrep.zip
  cp haplogrep.jar $JDIR/
fi

if [ ! -s $JDIR/mutserve.jar ] ; then
  wget https://github.com/seppinho/mutserve/releases/download/v1.3.4/mutserve-1.3.4.jar
  cp mutserve-1.3.4.jar $JDIR/mutserve.jar
  cd ../
fi

cd -
#################

if [ ! -s $RDIR/hs38DH.fa ] ; then
  cd $RDIR
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O hs38DH.fa
  $BDIR/samtools faidx hs38DH.fa
  cd -
fi

