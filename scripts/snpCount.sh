#!/usr/bin/env bash
set -e

if [ "$#" -ne 4 ]; then exit 0 ; fi

export IN=$1
D=$2
export M=$3 # source (mutect2...)
export T=$4 # thold
#export TL="env${T}pc_"
export TL=""

#######################################################

test -f $D/$M.00.concat.vcf

#catenate, merge, sount SNPs
cat $D/$M.00.concat.vcf | filterVcf.pl -p 0.$T > $D/$M.$T.concat.vcf
cat $D/$M.$T.concat.vcf | concat2cat.pl > $D/$M.$T.vcf
cat $D/$M.$T.concat.vcf | concat2merge.pl -in $IN | bedtools sort -header | tee $D/$M.$T.merge.vcf | vcf2sitesOnly.pl >  $D/$M.$T.merge.sitesOnly.vcf
cat $D/$M.$T.concat.vcf | snpCount.pl -in $IN | tee $D/$M.$T.tab | getSummaryN.pl > $D/$M.$T.summary

#######################################################
if [[ -z "${FNAME}" ]]; then exit 0;  fi
cat $D/$M.$T.concat.vcf | eval $FRULE > $D/$M.$T.$FNAME.concat.vcf
cat $D/$M.$T.$FNAME.concat.vcf | concat2merge.pl -in $IN | bedtools sort -header | tee $D/$M.$T.$FNAME.merge.vcf | vcf2sitesOnly.pl >  $D/$M.$T.$FNAME.merge.sitesOnly.vcf
cat $D/$M.$T.$FNAME.concat.vcf | snpCount.pl -in $IN | tee $D/$M.$T.$FNAME.tab | getSummaryN.pl > $D/$M.$T.$FNAME.summary

#cat $D/$M.$T.concat.vcf | uniqVcf.pl | posCount.pl  | sort -k1,1 -k2,2n | tee $D/$M.$T.count | awk 'NR == 1; NR > 1 {print $0 | "sort -k4,4nr -k2,2n"}' > $D/$M.$T.hcount
#cat $D/$M.$T.count |  awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6nr -k2,2n"}'  > $D/$M.$T.hScount
