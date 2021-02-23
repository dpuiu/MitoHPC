#!/bin/bash -eux

M=$1		# mutect2
T=$2		# 03
export F=$3	# NUMT
IN=${4:-in.txt}	# in.txt
I=0		# column

cat $M.$T.vcf  | grep -v -P ";$F" | cat2concat.pl | tee $M.$T.concat.no$F.vcf | snpCount.pl  -in $IN -i $I > $M.$T.no$F.tab

