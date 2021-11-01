#!/usr/bin/env bash
set -e

if [ "$#" -ne 2 ]; then exit 0 ; fi

export M=$1 # source (mutect2...)
export T=$2 # thold

#######################################################

test -f $HP_ODIR/$M.00.concat.vcf

#catenate, merge, sount SNPs
cat $HP_ODIR/$M.00.concat.vcf | filterVcf.pl -p 0.$T > $HP_ODIR/$M.$T.concat.vcf
cat $HP_ODIR/$M.$T.concat.vcf | concat2cat.pl > $HP_ODIR/$M.$T.vcf
cat $HP_ODIR/$M.$T.concat.vcf | concat2merge.pl -in $HP_IN | bedtools sort -header | tee $HP_ODIR/$M.$T.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$M.$T.merge.sitesOnly.vcf
cat $HP_ODIR/$M.$T.concat.vcf | snpCount.pl -in $HP_IN | tee $HP_ODIR/$M.$T.tab | getSummaryN.pl > $HP_ODIR/$M.$T.summary

#######################################################
if [[ -z "${FNAME}" ]]; then exit 0;  fi
cat $HP_ODIR/$M.$T.concat.vcf | eval $FRULE > $HP_ODIR/$M.$T.$FNAME.concat.vcf
cat $HP_ODIR/$M.$T.$FNAME.concat.vcf | concat2merge.pl -in $HP_IN | bedtools sort -header | tee $HP_ODIR/$M.$T.$FNAME.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$M.$T.$FNAME.merge.sitesOnly.vcf
cat $HP_ODIR/$M.$T.$FNAME.concat.vcf | snpCount.pl -in $HP_IN | tee $HP_ODIR/$M.$T.$FNAME.tab | getSummaryN.pl > $HP_ODIR/$M.$T.$FNAME.summary

#cat $HP_ODIR/$M.$T.concat.vcf | uniqVcf.pl | posCount.pl  | sort -k1,1 -k2,2n | tee $D/$M.$T.count | awk 'NR == 1; NR > 1 {print $0 | "sort -k4,4nr -k2,2n"}' > $HP_ODIR/$M.$T.hcount
#cat $HP_ODIR/$M.$T.count |  awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6nr -k2,2n"}'  > $HP_ODIR/$M.$T.hScount
