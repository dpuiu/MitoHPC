#!/usr/bin/env bash
set -ex

##########################################################

#Program that summarizes aligned read counts
#Input arguments

###########################################################
#get count stats

cut -f2 $HP_IN | sed -r "s|(.*)\.|\1\t|g" | cut -f1 | sed "s|$|.count|" | xargs cat | $HP_SDIR/uniq.pl | $HP_SDIR/getCN.pl > $HP_ODIR/count.tab
if [ $HP_I -lt 1 ] ; then exit 0 ; fi

###########################################################
# get 1st iteration stats

M=$HP_M

#cvg
cut -f3 $HP_IN | sed "s|$|.cvg.stat|" | xargs cat | uniq.pl -i 0  > $HP_ODIR/cvg.tab
cut -f3 $HP_IN | sed "s|$|.$M.00.vcf|" | xargs cat | bedtools sort -header  > $HP_ODIR/$M.00.concat.vcf
cat $HP_ODIR/$M.00.concat.vcf | grep "^#" > $HP_ODIR/$M.00.concat.vcf+
cat $HP_ODIR/$M.00.concat.vcf | grep -v "^#" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> $HP_ODIR/$M.00.concat.vcf+
mv  $HP_ODIR/$M.00.concat.vcf+  $HP_ODIR/$M.00.concat.vcf

#snv counts
snpCount.sh $M $HP_T1
snpCount.sh $M $HP_T2
snpCount.sh $M $HP_T3

if [[ ! -z "${HP_FNAME}" ]]; then
  cat $HP_ODIR/$M.00.concat.vcf | eval $HP_FRULE > $HP_ODIR/$M.$HP_FNAME.00.concat.vcf
  snpCount.sh $M.$HP_FNAME $HP_T1
  snpCount.sh $M.$HP_FNAME $HP_T2
  snpCount.sh $M.$HP_FNAME $HP_T3
fi

#cleanup
rm -f fastp.html fastp.json

#haplogroups
cut -f3 $HP_IN | sed "s|$|.$M.haplogroup|" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
  perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | sed "s|\.MT||" | \
  tee $HP_ODIR/$M.haplogroup.tab | \
  perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $HP_ODIR/$M.haplogroup1.tab

#haplocheck
cut -f3 $HP_IN | sed "s|$|.$M.haplocheck|" | xargs cat  | uniq.pl | sed 's|^"Sample"|"Run"|' | sed 's|"||g' | sed 's| ||g' > $HP_ODIR/$M.haplocheck.tab

#fasta
cut -f3 $HP_IN | sed "s|$|.$M.fa|"        | xargs cat > $HP_ODIR/$M.fa
samtools faidx  $HP_ODIR/$M.fa

cut -f3 $HP_IN | sed "s|$|.count|" | xargs cat | $HP_SDIR/uniq.pl > $HP_ODIR/count1.tab

##########################################################
# get 2nd iteration stats

if [ $HP_I -lt 2 ] ; then exit 0 ; fi
if [ $HP_M != "mutect2" ] ; then exit 0 ; fi

MM=$M.$M

#cvg
cut -f3 $HP_IN | sed "s|$|$M.cvg.stat|" | xargs cat | uniq.pl -i 0  > $HP_ODIR/$M.cvg.tab
cut -f3 $HP_IN | sed "s|$|.$MM.00.vcf|" | xargs cat | bedtools sort -header  > $HP_ODIR/$MM.00.concat.vcf
cat $HP_ODIR/$MM.00.concat.vcf | grep "^#" > $HP_ODIR/$MM.00.concat.vcf+
cat $HP_ODIR/$MM.00.concat.vcf | grep -v "^#" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> $HP_ODIR/$MM.00.concat.vcf+
mv  $HP_ODIR/$MM.00.concat.vcf+  $HP_ODIR/$MM.00.concat.vcf

#snv counts
snpCount.sh $MM $HP_T1
snpCount.sh $MM $HP_T2
snpCount.sh $MM $HP_T3

if [[ ! -z "${HP_FNAME}" ]]; then
  cat $HP_ODIR/$MM.00.concat.vcf | eval $HP_FRULE > $HP_ODIR/$MM.$HP_FNAME.00.concat.vcf
  snpCount.sh $MM.$HP_FNAME $HP_T1
  snpCount.sh $MM.$HP_FNAME $HP_T2
  snpCount.sh $MM.$HP_FNAME $HP_T3
fi

cut -f3 $HP_IN | sed "s|$|.$HP_M.count|" | xargs cat | $HP_SDIR/uniq.pl > $HP_ODIR/count2.$HP_M.tab
