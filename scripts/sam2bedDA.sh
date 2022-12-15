#!/bin/bash -eux

samtools view $1.bam chrM | sam2bedDA.pl -min 1500  -max 15069 | sort -u > $2.da.bed
cat $2.da.bed  |   cut -f1,3,9 | sort -V | uniq -c | perl -ane 'if($F[0]>=5) { print join "\t",@F[1,2,3,0]; print "\n"}' | bed2bed.pl > $2.da.count2_8.min5.bed
