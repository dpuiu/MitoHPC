#!/usr/bin/env bash

export N=$1
export S=$2
export E=$3
SE=$4

test -s $S.fa
test -s $S.fa.fai

samtools faidx $S.fa $S:1-$E | cat $S.fa - | grep -v ">" |  perl -ane 'BEGIN { print ">$ENV{N}\n" }  chomp; print ; END {print "\n"}'  > $SE.fa

bwa index $SE.fa -p $SE
samtools faidx $SE.fa

rm -f $SE.dict
java -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $SE.fa --OUTPUT $SE.dict
