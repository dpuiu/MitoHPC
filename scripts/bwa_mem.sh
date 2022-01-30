#!/usr/bin/env bash
set -ex

I=$1
O=$2

test -s $HP_RDIR/$HP_RNAME.fa

test -s $HP_RDIR/$HP_RNAME.amb
test -s $HP_RDIR/$HP_RNAME.ann
test -s $HP_RDIR/$HP_RNAME.bwt
test -s $HP_RDIR/$HP_RNAME.pac
test -s $HP_RDIR/$HP_RNAME.sa

if [ ! -s $O.bam ] ; then
  test -s ${I}_1.fq.gz
  test -s ${I}_2.fq.gz

  bwa mem $HP_RDIR/$HP_RNAME ${I}_[12].fq.gz -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view -bu | \
    samtools sort -m $HP_MM -@ $HP_P -o $O.bam
fi

if [ ! -s $I.bam.bai ] ; then
  samtools index $O.bam -@ $HP_P
fi

if [ ! -s $I.idxstats ] ; then
  samtools idxstats $O.bam > $O.idxstats
fi
