#!/bin/bash -eux

export SDIR=`dirname $0`        # script directory
source $SDIR/init_marcc.sh

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

test -f $JDIR/gatk.jar
test -f $JDIR/haplogrep.jar
test -f $JDIR/mutserve.jar
test -f $JDIR/picard.jar
test -f $RDIR/$H
test -f $RDIR/$R

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
java --version | head -1               >> checkInstall.log
 
echo "########################"  >> checkInstall.log
echo "JARS:"  >> checkInstall.log

ls -l $JDIR/picard.jar    | awk '{print $11}' >> checkInstall.log
ls -l $JDIR/gatk.jar      | awk '{print $11}' >> checkInstall.log
ls -l $JDIR/haplogrep.jar | awk '{print $11}' >> checkInstall.log
ls -l $JDIR/mutserve.jar  | awk '{print $11}' >> checkInstall.log

echo "########################"  >> checkInstall.log
echo "VARS:"  >> checkInstall.log

echo   "SDIR=" $SDIR >> checkInstall.log
echo   "JDIR=" $JDIR >> checkInstall.log
echo   "RDIR=" $RDIR >> checkInstall.log
echo   "H=" $H >> checkInstall.log
echo   "R=" $R >> checkInstall.log
echo   "M=" $M >> checkInstall.log

echo Success!
