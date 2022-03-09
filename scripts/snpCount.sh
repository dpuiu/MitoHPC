#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that generates SNV counts

##############################################################################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi

export M=$1 # source (mutect2...)
export T=$2 # thold

#######################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi
if [ -s $HP_ODIR/$M.00.concat.vcf ] ; then  cat $HP_ODIR/$M.00.concat.vcf | filterVcf.pl -p 0.$T > $HP_ODIR/$M.$T.concat.vcf ; fi

cat $HP_ODIR/$M.$T.concat.vcf | concat2cat.pl > $HP_ODIR/$M.$T.vcf
cat $HP_ODIR/$M.$T.concat.vcf | concat2merge.pl -in $HP_IN | bedtools sort -header | tee $HP_ODIR/$M.$T.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$M.$T.merge.sitesOnly.vcf
cat $HP_ODIR/$M.$T.concat.vcf | snpCount.pl -in $HP_IN | tee $HP_ODIR/$M.$T.tab | getSummaryN.pl > $HP_ODIR/$M.$T.summary

# get suspicious samples
if [ -f $HP_ODIR/$M.merge.bed ] ; then
  rm -f $HP_ODIR/$M.$T.suspicious.tab ; touch $HP_ODIR/$M.$T.suspicious.tab

  # low mtDNA_CN
  cat $HP_ODIR/count.tab | perl -ane 'print "$F[0]\tlow_CN\t$F[-1]\n" if($F[4] and $F[4]=~/^\d+/ and $F[4]<300/($ENV{T}+1)) ;' >> $HP_ODIR/$M.$T.suspicious.tab

  # mismatch haplogroup
  if [ -s $HP_ODIR/$M.haplogroup1.tab ] ; then
    sort $HP_ODIR/$M.haplogroup1.tab > $HP_ODIR/$M.haplogroup1.srt.tab
    cat $HP_ODIR/$M.$T.vcf | grep "AF=0" | grep "HG=" | perl -ane 'print "$F[-1]\t$1\n" if(/HG=(.+?);/)' | sort | uniq -c  | \
      perl -ane 'print "$F[1]\t$F[2]\t$F[0]\n" if($F[0]>1);'| sort  | join -  $HP_ODIR/$M.haplogroup1.srt.tab | \
      perl -ane 'print "$F[0]\tmismatch_HG\t$F[-1]\t$F[1]\t$F[2]\n" if($F[1] ne $F[-1]);' >> $HP_ODIR/$M.$T.suspicious.tab
    rm $HP_ODIR/$M.haplogroup1.srt.tab
  fi
  
  # multiple NUMT's
  cat $HP_ODIR/$M.$T.vcf | grep "AF=0" | grep "NUMT=" | perl -ane '/NUMT=(.+?);/ ; foreach(split /\|/,$1) { print "$F[-1]\t$_\n"}' | sort | uniq -c  | \
    perl -ane 'print "$F[1]\tmultiple_NUMTs\t$F[2]\t$F[0]\n" if($F[0]>1);' >> $HP_ODIR/$M.$T.suspicious.tab

  # haplockeck
  cat $HP_ODIR/$M.haplocheck.tab  | perl -ane 'print "$F[0]\thaplocheck_fail\t$F[2]\n" if($F[1] eq "YES" and $F[2]>$ENV{T}/100)' >> $HP_ODIR/$M.$T.suspicious.tab

  cut -f1 $HP_ODIR/$M.$T.suspicious.tab | sort -u > $HP_ODIR/$M.$T.suspicious.ids
fi

#######################################################
if [[ -z "${FNAME}" ]]; then exit 0;  fi
cat $HP_ODIR/$M.$T.concat.vcf | eval $FRULE > $HP_ODIR/$M.$T.$FNAME.concat.vcf
cat $HP_ODIR/$M.$T.$FNAME.concat.vcf | concat2merge.pl -in $HP_IN | bedtools sort -header | tee $HP_ODIR/$M.$T.$FNAME.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$M.$T.$FNAME.merge.sitesOnly.vcf
cat $HP_ODIR/$M.$T.$FNAME.concat.vcf | snpCount.pl -in $HP_IN | tee $HP_ODIR/$M.$T.$FNAME.tab | getSummaryN.pl > $HP_ODIR/$M.$T.$FNAME.summary

#cat $HP_ODIR/$M.$T.concat.vcf | uniqVcf.pl | posCount.pl  | sort -k1,1 -k2,2n | tee $D/$M.$T.count | awk 'NR == 1; NR > 1 {print $0 | "sort -k4,4nr -k2,2n"}' > $HP_ODIR/$M.$T.hcount
#cat $HP_ODIR/$M.$T.count |  awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6nr -k2,2n"}'  > $HP_ODIR/$M.$T.hScount
