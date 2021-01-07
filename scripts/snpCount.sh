#!/bin/bash -eux

##########################################################

#Program that summarizes aligned read counts
#Input arguments

D=$1 # out dir
M=$2 # snp caller: mutect2 or mutserve

find $D/ -name "*.$M.haplogroup" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' > $D/$M.haplogroup.tab
cat $D/$M.haplogroup.tab |  perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $D/$M.haplogroup1.tab

snpCount1.sh $D $M 03
snpCount1.sh $D $M 05
snpCount1.sh $D $M 10
snpCount1.sh $D $M max

join.pl  $D/$M.haplogroup.tab $D/$M.03.tab | join.pl - $D/$M.05.tab | join.pl - $D/$M.10.tab > $D/$M.tab
rm -f $D/$M.{03,05,10,max}.tab 
