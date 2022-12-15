#!/bin/bash -eux

samtools view $1.bam chrM | sam2bedSA.pl -min_score 140 -min_pos 125 -max_pos 16444  | sort -u > $2.sa.bed
cat $2.sa.bed  |   cut -f1,3,9 | sort -V | uniq -c | perl -ane 'if($F[0]>=5) { print join "\t",@F[1,2,3,0]; print "\n"}' | bed2bed.pl > $2.sa.count2_8.min5.bed
