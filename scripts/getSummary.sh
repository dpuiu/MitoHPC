#!/usr/bin/env bash
set -ex

##########################################################

#Program that summarizes aligned read counts
#Input arguments

###########################################################
#get count stats

if [ $HP_CN ]  && [ $HP_CN -ne 0 ] ; then awk '{print $3}' $HP_IN | sed "s|$|.all.count|" | xargs cat | uniq.pl | getCN.pl > $HP_ODIR/all.count.tab 
else                                      awk '{print $3}' $HP_IN | sed "s|$|.all.count|" | xargs cat | uniq.pl            > $HP_ODIR/all.count.tab
fi
if [ $HP_I -lt 1 ] ; then exit 0 ; fi

###########################################################
# get 1st iteration stats

M=$HP_M

#count,cvg
awk '{print $3}' $HP_IN | sed "s|$|.count|" | xargs cat | uniq.pl > $HP_ODIR/count.tab
awk '{print $3}' $HP_IN | sed "s|$|.cvg.stat|" | xargs cat | uniq.pl -i 0  > $HP_ODIR/cvg.tab
awk '{print $3}' $HP_IN | sed "s|$|.$M.00.vcf|" | xargs cat | bedtools sort -header  > $HP_ODIR/$M.00.concat.vcf
snpSort.sh $HP_ODIR/$M.00.concat
cat $HP_ODIR/$M.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $HP_ODIR/$M.00.AF.histo

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
if [ "$HP_O" == "Human" ] ; then
  awk '{print $3}' $HP_IN | sed "s|$|.$M.haplogroup|" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
    perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | sed "s|\.MT||" | \
    tee $HP_ODIR/$M.haplogroup.tab | \
    perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $HP_ODIR/$M.haplogroup1.tab

  #haplocheck
  awk '{print $3}' $HP_IN | sed "s|$|.$M.haplocheck|" | xargs cat  | uniq.pl | sed 's|^"Sample"|"Run"|' | sed 's|"||g' | sed 's| ||g' > $HP_ODIR/$M.haplocheck.tab
fi

#fasta
awk '{print $3}' $HP_IN | sed "s|$|.$M.fa|"        | xargs cat > $HP_ODIR/$M.fa
samtools faidx  $HP_ODIR/$M.fa

##########################################################
# get 2nd iteration stats

if [ $HP_I -lt 2 ] ; then exit 0 ; fi
if [ $HP_M != "mutect2" ] ; then exit 0 ; fi

MM=$M.$M
#count,cvg
awk '{print $3}' $HP_IN | sed "s|$|.$M.merge.bed|" | xargs cat > $HP_ODIR/$M.merge.bed
awk '{print $3}' $HP_IN | sed "s|$|.$M.count|" | xargs cat | uniq.pl > $HP_ODIR/$M.count.tab
awk '{print $3}' $HP_IN | sed "s|$|.$M.cvg.stat|" | xargs cat | uniq.pl -i 0  > $HP_ODIR/$M.cvg.tab
awk '{print $3}' $HP_IN | sed "s|$|.$MM.00.vcf|" | xargs cat | bedtools sort -header  > $HP_ODIR/$MM.00.concat.vcf
snpSort.sh $HP_ODIR/$MM.00.concat
cat $HP_ODIR/$MM.00.concat.vcf | grep -v "^#" | sed 's|:|\t|g'  | count.pl -i -1 -round 100| sort -n > $HP_ODIR/$MM.00.AF.histo

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
