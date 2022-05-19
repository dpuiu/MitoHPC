#!/usr/bin/env bash
set -eux

##############################################################################################################

# Program that downloads and install software prerequisites

##############################################################################################################

#. $HP_SDIR/init.sh

cd $HP_HDIR/
mkdir -p prerequisites/ $HP_BDIR/ $HP_JDIR/

cd prerequisites/
wget -N -c https://iweb.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 --no-check-certificate
wget -N -c https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
wget -N -c https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2
wget -N -c https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2
wget -N -c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
wget -N -c http://opengene.org/fastp/fastp
wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
wget -N -c https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip
wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip
wget -N -c https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip
wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13/mutserve.zip
wget -N -c https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz
