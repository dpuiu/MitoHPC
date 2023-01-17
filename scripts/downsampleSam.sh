#!/bin/bash -eux

IN=$1                                                                                   # input bam prefix(including path)
OUT=$2                                                                                  # output bam prefix(including path)
MCOUNT=${3:-10}                                                                         # max number of reads starting at a certain pos
R=${4:-"chrM chr1:629084-634672 chr17:22521208-22521639"}                               # regions: Ex chrM ...
MT=chrM

#Examples: 
#      downsampleSam.sh in out 
#      downsampleSam.sh in out 10
#      downsampleSam.sh in out 10 "chrM chr1:629084-634672 chr17:22521208-22521639"

test -s $IN.bam
test -s $IN.bam.bai
test -s $MT.fa.fai

################################

#downsample and index
samtools view -h $IN.bam $R -F 0x90C | filterSam.pl $R | downsampleSam.pl -max $MCOUNT | grep -v "^\@" | cut -f1 | sort | uniq -d  | samtools view -N /dev/stdin  $IN.bam -b  > $OUT.bam 
samtools index    $OUT.bam
samtools idxstats $OUT.bam > $OUT.idxstats

#compute cvg
samtools view -h $IN.bam  $MT | bedtools bamtobed | bedtools genomecov -g $MT.fa.fai -i /dev/stdin -d > $IN.cvg
samtools view -h $OUT.bam $MT | bedtools bamtobed | bedtools genomecov -g $MT.fa.fai -i /dev/stdin -d > $OUT.cvg
