#!/bin/bash -eux

D=$1
export M=$2 # source (mutect2...)
export T=$3 # thold

test -f $D/$M.haplogroup.tab
find $D/ -name "*.$M.$T.vcf" -not -name "*.$M.$M.$T.vcf"  | xargs cat | uniq.pl | bedtools sort -header > $D/$M.$T.vcf
if [ ! -s $D/$M.$T.vcf ] ; then exit 0 ; fi

annotateVcf.sh  $D/$M.$T.vcf

cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%H\n"} print;' > $D/$M.$T.H.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%h\n"} print;' > $D/$M.$T.h.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | grep "SNP;"    | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%S\n"} print;' > $D/$M.$T.S.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | grep "SNP;"    | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%s\n"} print;' > $D/$M.$T.s.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | grep "INDEL;"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%I\n"} print;' > $D/$M.$T.I.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | grep "INDEL;"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%i\n"} print;' > $D/$M.$T.i.tab

cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%Hp\n"} print;' > $D/$M.$T.Hp.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%hp\n"} print;' > $D/$M.$T.hp.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | grep "SNP;"    | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%Sp\n"} print;' > $D/$M.$T.Sp.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | grep "SNP;"    | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%sp\n"} print;' > $D/$M.$T.sp.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep -v ";AF=0" | grep "INDEL;"  | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%Ip\n"} print;' > $D/$M.$T.Ip.tab
cat $D/$M.$T.vcf | uniq2.pl -i 1 -j -1 | grep ";AF=0"    | grep "INDEL;"  | grep -v ";HP"  | count.pl -i 9  | sort | perl -ane 'BEGIN { print "Run\t$ENV{T}%ip\n"} print;' > $D/$M.$T.ip.tab

cut -f1  $D/$M.haplogroup.tab | \
  join.pl - $D/$M.$T.H.tab -empty 0  | join.pl - $D/$M.$T.h.tab -empty 0  | \
  join.pl - $D/$M.$T.S.tab -empty 0  | join.pl - $D/$M.$T.s.tab -empty 0  | \
  join.pl - $D/$M.$T.I.tab -empty 0  | join.pl - $D/$M.$T.i.tab -empty 0  | \
  join.pl - $D/$M.$T.Hp.tab -empty 0 | join.pl - $D/$M.$T.hp.tab -empty 0 | \
  join.pl - $D/$M.$T.Sp.tab -empty 0 | join.pl - $D/$M.$T.sp.tab -empty 0 | \
  join.pl - $D/$M.$T.Ip.tab -empty 0 | join.pl - $D/$M.$T.ip.tab -empty 0 > $D/$M.$T.tab

rm $D/$M.$T.?.tab $D/$M.$T.??.tab
