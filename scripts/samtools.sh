#!/usr/bin/env bash
set -ex

###############################################################

#Program that indexes and counts alignments in a BAM/CRAM file
#Input arguments:
#1: .bam/cram file; full path

################################################################

S=$1
test -s $2
N=`basename $2 .bam`
N=`basename $N .cram`

IDIR=`dirname $2`
I=$IDIR/$N
#PATH=$HP_SDIR:$HP_BDIR:$PATH

################################################################

#test BAM/CRAM file sorted
samtools view -H $2 | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

if [ ! -s $2.bai ] && [ ! -s $2.crai ]; then
  samtools index -@ $HP_P $2
fi

if [ ! -s $I.idxstats ] ; then
  samtools idxstats -@ $HP_P $2 > $I.idxstats
fi

if [ ! -s $I.count ] ; then
  cat $I.idxstats | idxstats2count.pl -sample $S -chrM $HP_MT >  $I.count
  samtools view -F 0x900 $2 $HP_NUMT -c | sed 's|^|NUMT\n|' | paste $I.count - > $I.count+ ; mv $I.count+ $I.count
fi
