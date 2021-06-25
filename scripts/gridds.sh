#!/bin/bash -eux

BDIR=/work-zfs/ssalzbe1/sw/packages/GRIDDS/
JDIR=/work-zfs/ssalzbe1/sw/packages/GRIDDS/
RDIR=~/Homo_sapiens_mito/RefSeq/
PATH=/work-zfs/ssalzbe1/sw/bin:$PATH

samtools view $1.bam chrM:200-16369 -b > $1.chrM_200_16369.bam
$BDIR/gridss -j $JDIR/gridss-2.12.0-gridss-jar-with-dependencies.jar -r $RDIR/chrM.fa -t 1 $1.chrM_200_16369.bam -o $1.gridds.vcf.gz -a $1.gridds.bam
rm $1.chrM_200_16369.bam

