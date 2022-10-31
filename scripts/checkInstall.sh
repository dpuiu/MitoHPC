#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that checks if all dependencies are installed

##############################################################################################################

#. $HP_SDIR/init.sh

#test executables and Java jars
which perl	        #usually available on Linux
which gcc	        #module load gcc
which java	        #module load java
which python		#module load python

which bwa	        #module load bwa
which samtools          #module load samtools
which bedtools          #module load bedtools
which fastp             #install from "https://github.com/OpenGene/fastp"
which samblaster        #install from "https://github.com/GregoryFaust/samblaster"
which bcftools
which tabix
which freebayes

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
#gridss

echo "########################"  >> checkInstall.log
echo "JAVA:"  >> checkInstall.log

#echo "JAVA_HOME=" $JAVA_HOME >> checkInstall.log
ls -l $HP_JDIR/gatk.jar       | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
ls -l $HP_JDIR/mutserve.jar   | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
ls -l $HP_JDIR/haplogrep.jar  | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
ls -l $HP_JDIR/haplocheck.jar | perl -ane 'print "$F[-1]\n";' >> checkInstall.log
#ls -l $HP_JDIR/gridss.jar

echo "########################"  >> checkInstall.log
echo "VARS:"  >> checkInstall.log

echo "HP_HDIR=" $HP_HDIR >> checkInstall.log
echo "HP_BDIR=" $HP_BDIR >> checkInstall.log
echo "HP_SDIR=" $HP_SDIR >> checkInstall.log
echo "HP_JDIR=" $HP_JDIR >> checkInstall.log
echo "HP_RDIR=" $HP_RDIR >> checkInstall.log

echo "########################"  >> checkInstall.log
echo "REFERENCES:"  >> checkInstall.log

test -s $HP_RDIR/$HP_RNAME.fa
test -s $HP_RDIR/$HP_MT.fa
test -s $HP_RDIR/$HP_MT+.fa
test -s $HP_RDIR/$HP_NUMT.fa

test -s $HP_RDIR/$HP_MT+.bwt
test -s $HP_RDIR/$HP_NUMT.bwt

test -s $HP_RDIR/$HP_MT.dict

echo Success!
