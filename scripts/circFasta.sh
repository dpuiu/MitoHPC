#!/usr/bin/env bash

test -s $S.fa
test -s $S.fa.fai

samtools faidx $S.fa $N:1-$E | cat $S.fa - | grep -v ">" |  perl -ane 'BEGIN { print ">$ENV{N}\n" }  chomp; print ; END {print "\n"}'  > $SE.fa

bwa index $SE.fa -p $SE
samtools faidx $SE.fa

rm -f $SE.dict
java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $SE.fa --OUTPUT $SE.dict
