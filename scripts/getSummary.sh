#!/usr/bin/env bash
set -e

#initialize ODIR (output dir) variable
if  [ "$#" -lt 1 ] ; then ODIR=$HP_ODIR ; else ODIR=$1 ; fi

##############################################################################################################

# Program that summarizes aligned read and SNV counts

##############################################################################################################

#get read count and mtDN-CN stats 
if [ $HP_CN ]  && [ $HP_CN -ne 0 ] || [ ! -s $ODIR/count.tab ] ; then
    cut -f2 $HP_IN | perl -ane 'print "$1.count\n" if(/(.+)\./);' | xargs cat | cut -f1,2,3,4 | uniq.pl | getCN.pl > $ODIR/count.tab
fi
if [ $HP_I -lt 1 ] ; then exit 0 ; fi

###########################################################
# get 1st iteration stats

S=$HP_M
V=$HP_V

#count,cvg
awk '{print $3}' $HP_IN | sed "s|$|.count|"     | xargs cat | uniq.pl > $ODIR/subsample.tab
awk '{print $3}' $HP_IN | sed "s|$|.cvg.stat|"  | xargs cat | uniq.pl -i 0  > $ODIR/cvg.tab
awk '{print $3}' $HP_IN | sed "s|$|.$S.00.vcf|" | xargs cat | uniq.pl | bedtools sort -header  > $ODIR/$S.00.concat.vcf
snpSort.sh $ODIR/$S.00.concat
cat $ODIR/$S.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $ODIR/$S.00.AF.histo
if [ $V ] ; then awk '{print $3}' $HP_IN | sed "s|$|.$V.00.vcf|" | xargs cat | uniq.pl | bedtools sort -header  > $ODIR/$V.00.concat.vcf ; fi

#haplogroups
if [ "$HP_O" == "Human" ] ; then
  awk '{print $3}' $HP_IN | sed "s|$|.$S.haplogroup|" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
    perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | sed "s|\.MT||" | \
    tee $ODIR/$S.haplogroup.tab | \
    perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $ODIR/$S.haplogroup1.tab

  #haplocheck
  awk '{print $3}' $HP_IN | sed "s|$|.$S.haplocheck|" | xargs cat  | uniq.pl | sed 's|^"Sample"|"Run"|' | sed 's|"||g' | sed 's| ||g' > $ODIR/$S.haplocheck.tab
fi

#fasta
awk '{print $3}' $HP_IN | sed "s|$|.$S.fa|"        | xargs cat > $ODIR/$S.fa
samtools faidx  $ODIR/$S.fa

#cvg
awk '{print $3}' $HP_IN | sed "s|$|.$S.merge.bed|" | xargs cat > $ODIR/$S.merge.bed

#snv counts
snpCount.sh $S $HP_T1
snpCount.sh $S $HP_T2
snpCount.sh $S $HP_T3

if [[ ! -z "${HP_FNAME}" ]]; then
  snpFilter.sh $S 00
  snpCount.sh $S.$HP_FNAME $HP_T1
  snpCount.sh $S.$HP_FNAME $HP_T2
  snpCount.sh $S.$HP_FNAME $HP_T3
fi
#cleanup
rm -f fastp.html fastp.json

##########################################################
# get 2nd iteration stats

if [ $HP_I -lt 2 ] ; then exit 0 ; fi
if [ $S == "mutserve" ] ; then exit 0 ; fi

SS=$S.$S
#count,cvg
awk '{print $3}' $HP_IN | sed "s|$|.$S.cvg.stat|" | xargs cat | uniq.pl -i 0  > $ODIR/$S.cvg.tab
awk '{print $3}' $HP_IN | sed "s|$|.$SS.00.vcf|"  | xargs cat | uniq.pl | bedtools sort -header  > $ODIR/$SS.00.concat.vcf
snpSort.sh $ODIR/$SS.00.concat
cat $ODIR/$SS.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $ODIR/$SS.00.AF.histo

#snv counts
snpCount.sh $SS $HP_T1
snpCount.sh $SS $HP_T2
snpCount.sh $SS $HP_T3

if [[ ! -z "${HP_FNAME}" ]]; then
  snpFilter.sh $SS 00
  snpCount.sh $SS.$HP_FNAME $HP_T1
  snpCount.sh $SS.$HP_FNAME $HP_T2
  snpCount.sh $SS.$HP_FNAME $HP_T3
fi


