#!/bin/bash -eux

D=$1
export M=$2 # source (mutect2...)
export T=$3 # thold

test -f $M.haplogroup.tab
find $D/ -name "*.$M.$T.vcf" | xargs cat | grep -v ^# | sort -k2,2n -k4,4 -k5,5 -k10,10  > $M.$T.vcf

cat $M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF" | grep "SNP;"    | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%S\n"} print;' > $M.$T.S.tab
cat $M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF"    | grep "SNP;"    | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%s\n"} print;' > $M.$T.s.tab

cat $M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF" | grep "INDEL;"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%I\n"} print;' > $M.$T.I.tab
cat $M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF"    | grep "INDEL;"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%i\n"} print;' > $M.$T.i.tab

cut -f1  $M.haplogroup.tab | \
  join.pl - $M.$T.S.tab -empty 0 | join.pl - $M.$T.s.tab -empty 0 | \
  join.pl - $M.$T.I.tab -empty 0 | join.pl - $M.$T.i.tab -empty 0 > $M.$T.tab

rm $M.$T.*.tab
