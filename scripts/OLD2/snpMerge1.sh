#!/bin/bash -eux

D=$1
export M=$2 # source (mutect2...)
export T=$3 # thold

###########################################################

find $D/ -name "*.$M.$T.vcf.gz" -not -name "*.$M.$M.$T.vcf.gz" | sort > $D/$M.$T.merge.txt
cat $D/$M.$T.merge.txt | xargs vcf-merge | fixmergeVcf.pl -in $D/$M.$T.merge.txt > $D/$M.$T.merge.vcf
bgzip -f $D/$M.$T.merge.vcf
tabix -f  $D/$M.$T.merge.vcf.gz

java -jar $JDIR/picard.jar  MakeSitesOnlyVcf --INPUT $D/$M.$T.merge.vcf.gz  --OUTPUT $D/$M.$T.merge.sitesOnly.vcf.gz 

