#!/usr/bin/env bash 

#######################################################################

#Program that downloads and install software prerequisites and hs38DH.fa

########################################################################


SDIR=`dirname $0`	 # script directory
source $SDIR/init.sh
mkdir -p $HDIR/prerequisites/ $BDIR/ $JDIR/

#############

cd $HDIR/prerequisites/
wget -N -c https://netactuate.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
wget -N -c https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
wget -N -c https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
wget -N -c https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
wget -N -c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
wget -N -c https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
wget -N -c https://github.com/nferraz/st/archive/v1.1.4.tar.gz -O st-1.1.4.tar.gz
wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
wget -N -c https://github.com/OpenGene/fastp/archive/v0.20.1.tar.gz -O fastp-0.20.1.tar.gz
#wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip 
wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip 
wget -N -c https://github.com/seppinho/mutserve/releases/download/v1.3.4/mutserve-1.3.4.jar 
wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc12/mutserve.zip 

cd $RDIR/ 
wget -N -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O hs38DH.fa

