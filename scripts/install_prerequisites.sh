#!/usr/bin/env bash
set -x

if [ -z $HP_SDIR ] ; then echo "Variable HP_SDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  exit 0 ; fi
if [ -z $HP_HDIR ] ; then echo "Variable HP_HDIR not defined. Make sure you followed the SETUP ENVIRONMENT instructions" ;  exit 0 ; fi

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
  wget -N -c --no-check-certificate https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
  if [ ! -s $HP_BDIR/bwa ] ; then
    tar -xjvf bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make ; cp bwa $HP_BDIR/
    cd -
  fi
fi

which samtools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2
  if [ ! -s $HP_BDIR/samtools ] ; then
    tar -xjvf samtools-1.16.tar.bz2
    cd samtools-1.16
    ./configure --prefix=$HP_HDIR/ # --without-curses --disable-bz2
    make ;  make install
    cd -
  fi
fi

which bcftools
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.16.tar.bz2
  if [ ! -s $HP_BDIR/bcftools ] ; then
    tar -xjvf  bcftools-1.16.tar.bz2
    cd bcftools-1.16
    ./configure --prefix=$HP_HDIR/ # --disable-bz2
    make ; make install
    cd -
  fi
fi

which htsfile
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget -N -c https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
  if [ ! -s $HP_BDIR/tabix ] ; then
    tar -xjvf htslib-1.16.tar.bz2
    cd htslib-1.16
    ./configure --prefix=$HP_HDIR/ # --disable-bz2
    make ; make install
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

  #wget -N -c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
  #cp bedtools.static.binary $HP_BDIR/bedtools
  #chmod a+x $HP_BDIR/bedtools
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

which gridss
if [[ $? != 0 || $# == 1 && $1 == "-f" ]] ; then
  wget https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz
  tar -xzvf gridss-2.13.2.tar.gz
  cp gridss $HP_BDIR/
  cp gridss-2.13.2-gridss-jar-with-dependencies.jar $HP_JDIR/gridss.jar
fi

if [ ! -s $HP_JDIR/gatk.jar ] ; then
  wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
  unzip -o gatk-4.3.0.0.zip
  cp gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar $HP_JDIR/gatk.jar
  cp gatk-4.3.0.0/gatk $HP_BDIR/
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
  wget -N -c https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip
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

if [ ! -s $HP_RDIR/$HP_MTC.fa ] ; then
  circFasta.sh $HP_MT $HP_E $HP_MTC
fi

if [ ! -s $HP_RDIR/$HP_MTR.fa ] ; then
  rotateFasta.sh $HP_MT $HP_E $HP_MTR
fi

