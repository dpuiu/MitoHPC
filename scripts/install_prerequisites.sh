#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that downloads and installs software prerequisites and genome reference
#  -f : opt; force reinstall
##############################################################################################################


#. $HP_SDIR/init.sh
cd $HP_HDIR
mkdir -p prerequisites/ $HP_BDIR/ $HP_JDIR/ $HP_RDIR/
cd prerequisites/

which bwa
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c --no-check-certificate https://iweb.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
  if [ ! -s $HP_BDIR/bwa ] ; then
    tar -xjvf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make ; cp bwa $HP_BDIR/
    cd -
  fi
fi

which samtools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
  if [ ! -s $HP_BDIR/samtools ] ; then
    tar -xjvf samtools-1.14.tar.bz2
    cd samtools-1.14
    ./configure --prefix=$HP_HDIR/ ; make ;  make install
    cd -
  fi
fi

which bcftools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
  if [ ! -s $HP_BDIR/bcftools ] ; then
    tar -xjvf  bcftools-1.14.tar.bz2
    cd bcftools-1.14
    ./configure --prefix=$HP_HDIR/ ; make ; make install
    cd -
  fi
fi

which htsfile
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
  if [ ! -s $HP_BDIR/tabix ] ; then
    tar -xjvf htslib-1.14.tar.bz2
    cd htslib-1.14
    ./configure --prefix=$HP_HDIR/ ; make ; make install
    cd -
  fi
fi

which samblaster
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
  if [ ! -s $HP_BDIR/samblaster ] ; then
    tar -xzvf samblaster-v.0.1.26.tar.gz
    cd samblaster-v.0.1.26
    make ; cp samblaster $HP_BDIR/
    cd -
  fi
fi

which bedtools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
  if [ ! -s $HP_BDIR/bedtools ] ; then
    tar -xzvf bedtools-2.30.0.tar.gz
    cd bedtools2/
    make install prefix=$HP_HDIR/
    cd -
  fi
fi

which fastp
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c http://opengene.org/fastp/fastp
  cp fastp $HP_BDIR/
  chmod a+x $HP_BDIR/fastp
fi

which freebayes
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
  gunzip freebayes-1.3.6-linux-amd64-static.gz  -c >  $HP_BDIR/freebayes
  chmod a+x $HP_BDIR//freebayes
fi

if [ ! -s $HP_JDIR/gatk.jar ] ; then
  wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.2.4.0/gatk-4.2.4.0.zip
  unzip -o gatk-4.2.4.0.zip
  cp gatk-4.2.4.0/gatk-package-4.2.4.0-local.jar $HP_JDIR/gatk.jar
  cp gatk-4.2.4.0/gatk $HP_BDIR/
fi

if [ ! -s $HP_JDIR/haplogrep.jar ] ; then
  wget -N -c https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip
  unzip -o haplogrep.zip
  cp haplogrep.jar $HP_JDIR/
fi

if [ ! -s $HP_JDIR/haplocheck.jar ] ; then
  wget -N -c https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip
  unzip -o haplocheck.zip
  cp haplocheck.jar $HP_JDIR/
fi

if [ ! -s $HP_JDIR/mutserve.jar ] ; then
  wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc12/mutserve.zip
  unzip -o mutserve.zip
  cp mutserve.jar $HP_JDIR
fi

if [ ! -s $HP_RDIR/$HP_RNAME.fa ] ; then
  wget -qO- $HP_RURL | zcat -f > $HP_RDIR/$HP_RNAME.fa
  samtools faidx $HP_RDIR/$HP_RNAME.fa
fi

if [ ! -s $HP_RDIR/$HP_MT.fa ] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RMT > $HP_RDIR/$HP_MT.fa
  samtools faidx $HP_RDIR/$HP_MT.fa
  java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $HP_RDIR/$HP_MT.fa --OUTPUT $HP_RDIR/$HP_MT.dict
fi

if [ ! -s $HP_RDIR/$HP_NUMT.fa ] ; then
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RNUMT > $HP_RDIR/$HP_NUMT.fa
  bwa index $HP_RDIR/$HP_NUMT.fa -p $HP_RDIR/$HP_NUMT
fi

if [ ! -s $HP_RDIR/${HP_MT}+.fa ] ; then
  cat $HP_RDIR/$HP_MT.fa > $HP_RDIR/${HP_MT}+.fa
  samtools faidx $HP_RDIR/$HP_RNAME.fa $HP_RMT:1-$HP_E | grep -v ">" >> $HP_RDIR/${HP_MT}+.fa
  #cat $HP_RDIR/$HP_NUMT.fa >> $HP_RDIR/${HP_MT}+.fa
  bwa index $HP_RDIR/${HP_MT}+.fa -p $HP_RDIR/${HP_MT}+
fi

