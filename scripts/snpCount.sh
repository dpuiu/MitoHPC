#!/bin/bash -eux

D=$1 # out dir
M=$2 # source(mutect2...)

if [ ! -s $M.haplogroup.tab ] ; then
  find $D/ -name "*.$M.haplogroup"  | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort | perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' > $M.haplogroup.tab
fi

snpCount1.sh $D $M 01
snpCount1.sh $D $M 03
snpCount1.sh $D $M 05
snpCount1.sh $D $M 10

join.pl  $M.haplogroup.tab $M.01.tab | join.pl - $M.03.tab | join.pl - $M.05.tab | join.pl - $M.10.tab | column -t > $M.tab
rm $M.??.tab
