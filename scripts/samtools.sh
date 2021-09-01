#!/usr/bin/env bash
set -e

###############################################################

#Program that indexes and counts alignments in a BAM/CRAM file
#Input arguments:
I=$1
#MT environmental variable must be set

################################################################

export N=`basename $I .bam`
export N=`basename $N .cram`
export D=`dirname $I`
export PATH=$HP_SDIR:$HP_BDIR:$PATH

P=1

test -f $I
################################################################

#test BAM/CRAM file sorted
samtools view -H $I | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

if [ ! -s $I.bai ] && [ ! -s $I.crai ]; then
  samtools index -@ $P $I
fi

if [ ! -s $D/$N.idxstats ] ; then
  samtools idxstats -@ $P $I > $D/$N.idxstats
fi

if [ ! -s $D/$N.count ] ; then
  cat $D/$N.idxstats | idxstats2count.pl -sample $N -chrM $HP_MT >  $D/$N.count
  samtools view -F 0x900 $I $HP_NUMT -c | sed 's|^|NUMT\n|' | paste $D/$N.count - > $D/$N.count+
  mv $D/$N.count+ $D/$N.count
fi
