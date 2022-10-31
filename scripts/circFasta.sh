#!/usr/bin/bash -eux

export N=$1
export S=$2
export E=$3

SE=$S.c$E

test -s $S.fa
test -s $S.fa.fai

samtools faidx $S.fa $S:1-$E | cat $S.fa - | grep -v ">" |  perl -ane 'BEGIN { print ">$ENV{N}\n" }  chomp; print ; END {print "\n"}'  > $SE.fa

bwa index $SE.fa -p $SE
samtools faidx $SE.fa


