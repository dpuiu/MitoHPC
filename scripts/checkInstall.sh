#!/usr/bin/env bash
set -e

#########################################################

#Program that checks if all dependencies are installed
#No input arguments

##########################################################

export HP_SDIR=`dirname $0`        # script directory
export HP_HDIR=`readlink -f $HP_SDIR/..`
export HP_LDIR=$HP_HDIR/lib/perl5
export PERLLIB=$HP_LDIR:$PERLLIB
export PERL5LIB=$HP_LDIR:$PERL5LIB
export HP_JDIR=$HP_HDIR/java
export HP_BDIR=$HP_HDIR/bin
export HP_RDIR=$HP_HDIR/RefSeq

export PATH=$HP_SDIR:$HP_BDIR:$PATH
export PERLLIB=$HP_LDIR/:$PERLLIB
export PERL5LIB=$HP_LDIR/:$PERL5LIB

#test executables and Java jars
which perl	        #usually available on Linux
which gcc	        #MARCC:"module load gcc"
which java	        #MARCC:"module load java

which bwa	        #MARCC:"module load bwa"  or install from "https://sourceforge.net/projects/bio-bwa/files/"
which samtools          #install from "http://www.htslib.org/download/"
which bedtools          #MARCC:"module load bedtools"
which fastp             #install from "https://github.com/OpenGene/fastp"
which samblaster        #install from "https://github.com/GregoryFaust/samblaster"
which circSam.pl        #available under scripts
which filterVcf.pl      #available under scripts
which st                #simple statistics: install from "https://github.com/nferraz/st"
perl -e 'use App::St'
which bcftools
which tabix
which vcftools

test -f $HP_JDIR/gatk.jar
test -f $HP_JDIR/mutserve.jar
test -f $HP_JDIR/haplogrep.jar
test -f $HP_JDIR/haplocheck.jar
######################################################

echo "########################"  > checkInstall.log
echo "DATE:"  >> checkInstall.log
date >> checkInstall.log

echo "########################"  >> checkInstall.log
echo "EXECUTABLES:"  >> checkInstall.log
samtools --version | head -1 >> checkInstall.log
bcftools --version | head -1 >> checkInstall.log
bedtools --version | head -1 >> checkInstall.log
fastp --version 2>&1 | head -1  >> checkInstall.log
echo -n "bwa "  >> checkInstall.log; bwa 2>&1 | grep Version -m 1  >> checkInstall.log
samblaster --version 2>&1 | head -1 >> checkInstall.log
echo -n "tabix " >> checkInstall.log ; tabix 2>&1 |  grep -v ^$ | head -1 >> checkInstall.log
perl --version | grep -v ^$ | head -1  >> checkInstall.log
java -version | head -1               >> checkInstall.log
 
echo "########################"  >> checkInstall.log
echo "JAVA_HOME=" $JAVA_HOME >> checkInstall.log
echo "JARS:"  >> checkInstall.log

ls -l $HP_JDIR/gatk.jar      | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
ls -l $HP_JDIR/haplogrep.jar | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
ls -l $HP_JDIR/mutserve.jar  | perl -ane 'print "$F[-1]\n";' >> checkInstall.log

echo "########################"  >> checkInstall.log
echo "VARS:"  >> checkInstall.log

echo "HP_HDIR=" $HP_HDIR >> checkInstall.log
echo "HP_BDIR=" $HP_BDIR >> checkInstall.log
echo "HP_SDIR=" $HP_SDIR >> checkInstall.log
echo "HP_JDIR=" $HP_JDIR >> checkInstall.log
echo "HP_RDIR=" $HP_RDIR >> checkInstall.log

echo "########################"  >> checkInstall.log
echo "REFERENCE SEQUENCES:"  >> checkInstall.log

#test -s $HP_RDIR/$HP_HG
#test -s $HP_RDIR/$HP_MT.fa
#test -s $HP_RDIR/$HP_R

echo Success!
