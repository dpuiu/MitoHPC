#!/usr/bin/env bash 

##########################################################

#Program that summarizes aligned read counts
#Input arguments

D=$1 			# out dir
M=$2 			# snp caller: mutect2 or mutserve
T1=${3:-03}              # Heteroplasmy tholds
T2=${4:-05}
T3=${5:-10}

find $D/ -name "*.$M.haplogroup" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' > $D/$M.haplogroup.tab
cat $D/$M.haplogroup.tab |  perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $D/$M.haplogroup1.tab

snpCount1.sh $D $M $T1
snpCount1.sh $D $M $T2
snpCount1.sh $D $M $T3
#snpCount1.sh $D $M max

join.pl  $D/$M.haplogroup.tab $D/$M.$T1.tab | join.pl - $D/$M.$T2.tab | join.pl - $D/$M.$T3.tab > $D/$M.tab
rm -f $D/$M.{$T1,$T2,$T3,max}.tab 
