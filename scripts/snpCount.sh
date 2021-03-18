#!/usr/bin/bash -eux


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
cat $D/$M.$T.concat.vcf | concat2merge.pl -in $IN | bedtools sort -header | tee $D/$M.$T.merge.vcf | vcf2sitesOnly.pl >  $D/$M.$T.merge.sitesOnly.vcf
cat $D/$M.$T.concat.vcf | snpCount.pl -in $IN | tee $D/$M.$T.tab | getSummaryN.pl > $D/$M.$T.summary

#######################################################
if [[ -z "${FNAME}" ]]; then exit 0;  fi

cat $D/$M.$T.concat.vcf | $FRULE > $D/$M.$T.$FNAME.concat.vcf
cat $D/$M.$T.$FNAME.concat.vcf | concat2merge.pl -in $IN | bedtools sort -header | tee $D/$M.$T.$FNAME.merge.vcf | vcf2sitesOnly.pl >  $D/$M.$T.$FNAME.merge.sitesOnly.vcf
cat $D/$M.$T.$FNAME.concat.vcf | snpCount.pl -in $IN | tee $D/$M.$T.$FNAME.tab | getSummaryN.pl > $D/$M.$T.$FNAME.summary

