#!/bin/bash -eux

IN=$1                                                                                   # input bam prefix(including path)
OUT=$2                                                                                  # output bam prefix(including path)
MCOUNT=${3:-10}                                                                         # max number of reads starting at a certain pos
R=${4:-"chrM chr1:629084-634672 chr17:22521208-22521639"}                               # MT+NUMT regions (MT 1st) 
                                                                                        # chr1:76971223-76971280 & chr17:22521208-22521639 align to chrM 5'/3'
                                                                                        # chr1:629084-634672 : 5Kb NUMT

#Examples: 
#      downsampleSam.sh in out 
#      downsampleSam.sh in out 10
#      downsampleSam.sh in out 10 chrM
#      downsampleSam.sh in out 10 "chrM chr1:629084-634672 chr1:76971223-76971280 chr17:22521208-22521639"

test -s $IN.bam
#test -s $IN.bam.bai

################################

#downsample and index
if [ ! -f  $OUT.bam ] ; then
  samtools view -h $IN.bam $R -F 0x90C | filterSam.pl $R | \
    samtools view -b | samtools sort -n | samtools fixmate /dev/stdin /dev/stdout | samtools view -h | \
    downsampleSam.pl -max $MCOUNT | grep -v "^\@" | cut -f1 | sort | uniq -d  | samtools view -N /dev/stdin  $IN.bam -b > $OUT.bam
  samtools index    $OUT.bam
  samtools idxstats $OUT.bam > $OUT.idxstats
fi
