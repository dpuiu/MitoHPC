#!/usr/bin/bash -eux

export S=$1
export E=$2
SE=$3

test -s $S.fa
test -s $S.fa.fai

samtools faidx $S.fa $S:1-$E | cat $S.fa - | grep -v ">" |  perl -ane 'BEGIN { print ">$ENV{S}\n" }  chomp; print ; END {print "\n"}'  > $SE.fa

bwa index $SE.fa -p $SE
samtools faidx $SE.fa


