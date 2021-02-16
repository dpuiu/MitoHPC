#!/bin/bash -eux

D=$1
export M=$2 # source (mutect2...)
export T=$3 # thold
export TL="env${T}pc_"

test -f $D/$M.haplogroup.tab

#######################################################
#catenate

find $D/ -name "*.$M.$T.vcf.gz" -not -name "*.$M.$M.$T.vcf.gz" | sort > $D/$M.$T.merge.txt

cat $D/$M.$T.merge.txt | xargs vcf-concat | bedtools sort -header >  $D/$M.$T.concat.vcf
annotateVcf.sh  $D/$M.$T.concat.vcf
bgzip -f $D/$M.$T.concat.vcf

#########################################################
#count

zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    |                 perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}H\n"} print;' > $D/$M.$T.H.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" |                 perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}h\n"} print;' > $D/$M.$T.h.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    | grep -v ";INDEL"   | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}S\n"} print;' > $D/$M.$T.S.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" | grep -v ";INDEL"   | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}s\n"} print;' > $D/$M.$T.s.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    | grep ";INDEL" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}I\n"} print;' > $D/$M.$T.I.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" | grep ";INDEL" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}i\n"} print;' > $D/$M.$T.i.tab

zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    | grep -v ";HP" |                 perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}Hp\n"} print;' > $D/$M.$T.Hp.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" | grep -v ";HP" |                 perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}hp\n"} print;' > $D/$M.$T.hp.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    | grep -v ";INDEL"   | grep -v ";HP" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}Sp\n"} print;' > $D/$M.$T.Sp.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" | grep -v ";INDEL"   | grep -v ";HP" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}sp\n"} print;' > $D/$M.$T.sp.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep ":1$"    | grep ";INDEL" | grep -v ";HP" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}Ip\n"} print;' > $D/$M.$T.Ip.tab
zcat $D/$M.$T.concat.vcf.gz | uniqSnp.pl | grep -v "^#" | grep -v ":1$" | grep ";INDEL" | grep -v ";HP" | perl -ane 'print "$1\n" if(/SM=(.+?)[;\s]/);' | count.pl | sort | perl -ane 'BEGIN { print "Run\t$ENV{TL}ip\n"} print;' > $D/$M.$T.ip.tab

cut -f1  $D/$M.haplogroup.tab | \
  join.pl - $D/$M.$T.H.tab -empty 0  | join.pl - $D/$M.$T.h.tab -empty 0  | \
  join.pl - $D/$M.$T.S.tab -empty 0  | join.pl - $D/$M.$T.s.tab -empty 0  | \
  join.pl - $D/$M.$T.I.tab -empty 0  | join.pl - $D/$M.$T.i.tab -empty 0  | \
  join.pl - $D/$M.$T.Hp.tab -empty 0 | join.pl - $D/$M.$T.hp.tab -empty 0 | \
  join.pl - $D/$M.$T.Sp.tab -empty 0 | join.pl - $D/$M.$T.sp.tab -empty 0 | \
  join.pl - $D/$M.$T.Ip.tab -empty 0 | join.pl - $D/$M.$T.ip.tab -empty 0 > $D/$M.$T.tab

rm -f $D/$M.$T.?.tab $D/$M.$T.??.tab $D/$M.$T.concat.vcf.gz.tbi

#################################################################
#merge

cat $D/$M.$T.merge.txt | xargs vcf-merge | fixmergeVcf.pl -in $D/$M.$T.merge.txt > $D/$M.$T.merge.vcf
bgzip -f $D/$M.$T.merge.vcf
tabix -f  $D/$M.$T.merge.vcf.gz

java -jar $JDIR/picard.jar  MakeSitesOnlyVcf --INPUT $D/$M.$T.merge.vcf.gz  --OUTPUT $D/$M.$T.merge.sitesOnly.vcf.gz 
tabix -f $D/$M.$T.merge.sitesOnly.vcf.gz

