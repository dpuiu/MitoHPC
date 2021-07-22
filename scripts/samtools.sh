#!/usr/bin/env bash
set -ex

###############################################################

#Program that indexes and counts alignments in a BAM/CRAM file
#Input arguments:
I=$1
#MT environmental variable must be set

################################################################

export N=`basename $I .bam`
export N=`basename $N .cram`
export D=`dirname $I`
P=1

test -f $I
################################################################

#test BAM/CRAM file sorted 
samtools view -H $I | grep -m 1 -P "^@HD.+coordinate$"

if [ ! -s $I.bai ] && [ ! -s $I.crai ]; then
  samtools index -@ $P $I
fi

if [ ! -s $D/$N.idxstats ] ; then
  samtools idxstats -@ $P $I > $D/$N.idxstats
fi

if [ ! -s $D/$N.count ] ; then
  cat $D/$N.idxstats | idxstats2count.pl -sample $N -chrM $HP_MT >  $D/$N.count
fi

