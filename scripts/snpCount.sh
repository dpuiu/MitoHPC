#!/usr/bin/env bash
set -ex

##############################################################################################################

# Program that generates SNV counts

##############################################################################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi

export S=$1 # source (mutect2...)
export T=$2 # thold

if [ "$T" -eq 0 ]; then exit 0 ; fi

#######################################################

# 1st iteration;  get suspicious samples
if [ -f $HP_ODIR/$S.merge.bed ] ; then
  rm -f $HP_ODIR/$S.$T.suspicious.tab ; touch $HP_ODIR/$S.$T.suspicious.tab

  # low cvg (2022/06/01,10)
  cat $HP_ODIR/cvg.tab | perl -ane 'print "$F[0]\tmin_cvg_less_100\n"  if($F[7] and $F[7]=~/^\d+/ and $F[2]<100);'  >> $HP_ODIR/$S.$T.suspicious.tab
  cat $HP_ODIR/cvg.tab | perl -ane 'print "$F[0]\tmean_cvg_less_500\n"  if($F[7] and $F[7]=~/^\d+/ and $F[7]<500);'  >> $HP_ODIR/$S.$T.suspicious.tab

  # low mtDNA_CN
  cat $HP_ODIR/count.tab | perl -ane 'print "$F[0]\tlow_CN\t$F[-1]\n" if($F[4] and $F[4]=~/^\d+/ and $F[4]<300/($ENV{T}+1)) ;' >> $HP_ODIR/$S.$T.suspicious.tab

  # mismatch haplogroup
  if [ -s $HP_ODIR/$S.haplogroup1.tab ] ; then
    sort $HP_ODIR/$S.haplogroup1.tab > $HP_ODIR/$S.haplogroup1.srt.tab
    cat $HP_ODIR/$S.$T.vcf | grep "AF=0.0" | grep "HG=" | perl -ane 'print "$F[-1]\t$1\n" if(/HG=(.+?);/)' | sort | uniq -c  | \
      perl -ane 'print "$F[1]\t$F[2]\t$F[0]\n" if($F[0]>1);'| sort  | join -  $HP_ODIR/$S.haplogroup1.srt.tab | \
      perl -ane 'print "$F[0]\tmismatch_HG\t$F[-1]\t$F[1]\t$F[2]\n" if($F[1] ne $F[-1]);' >> $HP_ODIR/$S.$T.suspicious.tab
    rm $HP_ODIR/$S.haplogroup1.srt.tab
  fi
 
  # multiple low frequency NUMT's
  cat $HP_ODIR/$S.$T.vcf | grep "AF=0.0" | grep "NUMT=" | perl -ane '/NUMT=(.+?);/ ; foreach(split /\|/,$1) { print "$F[-1]\t$_\n"}' | sort | uniq -c  | \
    perl -ane 'print "$F[1]\tmultiple_NUMTs\t$F[2]\t$F[0]\n" if($F[0]>1);' >> $HP_ODIR/$S.$T.suspicious.tab

  # haplockeck
  cat $HP_ODIR/$S.haplocheck.tab  | perl -ane 'print "$F[0]\thaplocheck_fail\t$F[2]\n" if($F[1] eq "YES" and $F[2]>$ENV{T}/100)' >> $HP_ODIR/$S.$T.suspicious.tab

  cut -f1 $HP_ODIR/$S.$T.suspicious.tab | sort -u > $HP_ODIR/$S.$T.suspicious.ids
fi

########################################################

#Sep 22nd 2023; added "| eval $HP_FRULE"
#if [ -s $HP_ODIR/$S.00.concat.vcf ] ; then cat $HP_ODIR/$S.00.concat.vcf | filterVcf.pl -p 0.$T -suspicious $HP_ODIR/$HP_M.$T.suspicious.ids > $HP_ODIR/$S.$T.concat.vcf ; fi
if [ -s $HP_ODIR/$S.00.concat.vcf ] ; then cat $HP_ODIR/$S.00.concat.vcf | filterVcf.pl -p 0.$T -suspicious $HP_ODIR/$HP_M.$T.suspicious.ids | eval $HP_FRULE > $HP_ODIR/$S.$T.concat.vcf ; fi

test -s $HP_ODIR/$S.$T.concat.vcf

cat $HP_ODIR/$S.$T.concat.vcf | concat2cat.pl  > $HP_ODIR/$S.$T.vcf

# Jan 25 2023
#cat $HP_ODIR/$S.$T.concat.vcf | concat2merge.pl -in $HP_IN -suspicious $HP_ODIR/$HP_M.$T.suspicious.ids | bedtools sort -header | tee $HP_ODIR/$S.$T.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$S.$T.merge.sitesOnly.vcf 
#cat $HP_ODIR/$S.$T.concat.vcf | snpCount.pl     -in $HP_IN -suspicious $HP_ODIR/$HP_M.$T.suspicious.ids | tee $HP_ODIR/$S.$T.tab | getSummaryN.pl > $HP_ODIR/$S.$T.summary
#cat $HP_ODIR/$S.$T.concat.vcf | concat2pos.pl   -in $HP_IN -suspicious $HP_ODIR/$HP_M.$T.suspicious.ids | sort -k2,2n -k4,4 -k5,5  > $HP_ODIR/$S.$T.pos


# March 7th 2023
#cat $HP_ODIR/$S.$T.concat.vcf | concat2merge.pl -in $HP_IN  | bedtools sort -header | tee $HP_ODIR/$S.$T.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$S.$T.merge.sitesOnly.vcf 
cat $HP_ODIR/$S.$T.concat.vcf | concat2merge.pl -in $HP_IN  | tee $HP_ODIR/$S.$T.merge.vcf | vcf2sitesOnly.pl >  $HP_ODIR/$S.$T.merge.sitesOnly.vcf
cat $HP_ODIR/$S.$T.concat.vcf | snpCount.pl     -in $HP_IN  | tee $HP_ODIR/$S.$T.tab | getSummaryN.pl > $HP_ODIR/$S.$T.summary
cat $HP_ODIR/$S.$T.concat.vcf | concat2pos.pl   -in $HP_IN  | sort -k2,2n -k4,4 -k5,5  > $HP_ODIR/$S.$T.pos


