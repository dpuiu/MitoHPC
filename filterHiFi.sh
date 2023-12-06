#!/usr/bin/env bash
set -ex

#########################################################################################################################################
# Program that runs the heteroplasmy pipeline on a single HiFi sample

# Input arguments
#  1: sample names
#  2: BAM/CRAM alignment file; full path
#  3: output prefix; full path
#########################################################################################################################################

#set variables

export S=$1             # sample name
IDIR=`dirname "$2"`     # dir path
I=${2%.*}               # bam/cram file name prefix
X=${2##*.}              #                    extension
O=$3                    # output prefix
OS=$O.$HP_M             # output prefix + snv_caller
OSS=$OS.$HP_M
RG="@RG\tID:$S\tSM:$S\tPL:HiFi"
export PC="0.95"

ODIR=`dirname "$3"`; mkdir -p $ODIR

#########################################################################################################################################
# test if count and VCF output files exist; exit if they do

if [ $HP_I -lt 1 ] && [ -s $O.count ]    ; then exit 0 ; fi
if [ $HP_I -eq 1 ] && [ -s $OS.00.vcf  ] ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $OSS.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################
# test alignment file exists, is indexed and sorted by coordinates
test -s $2
test -s $2.bai || test -s $2.crai || test -s $I.bai || test -s $I.crai

if [ ! -s $O.fa ] ; then
  if [ "$X" == "cram" ] ; then T="-T $HP_RFILE.fa" ; else T="" ; fi
  MTCOUNT=`samtools view -F 0x900 $2 $HP_RMT $T -c`
  COUNT=`samtools view -F 0x900 $2 $T -c`
  echo -e "$S\t$COUNT\t$COUNT\$MTCOUNT" > $O.count

  if [ $MTCOUNT -lt 1 ]; then echo "ERROR: There are no MT reads in $2; plese remove it from $HP_IN" ; exit 1 ; fi
  if [ $HP_L ]; then R=`tail -1 $O.count | perl -ane '$R=$ENV{HP_L}/($F[-1]+1);  if($R<1) { print "-s $R"} else {print ""} '` ; else R="" ; fi

  samtools view -F 0x900 $R $2 $HP_RMT $T | perl -ane 'print "$F[0]\t",length($F[9]),"\n";' | sort > $O.len
  samtools view -F 0x900 $R $2 $HP_RMT $T -b | samtools fasta > $O.fa
fi

if [ ! -s $O.bam ] ; then
  cat $O.fa | minimap2  $HP_RDIR/$HP_MT.fa /dev/stdin -R $RG -a | samtools view -b | samtools sort | samtools view -h > $O.sam
  bedtools bamtobed -i $O.sam | cut -f 1-4 | bed2bed.pl  | count.pl -i 3 -j 4  |  sort | join $O.len -  -a 1 --nocheck-order | \
    perl -ane 'print "$F[0]\t$F[1]\t$F[2]\t",$F[2]/$F[1],"\n"' | perl -ane 'print  if($F[-1]>$ENV{PC});' | cut -f1 | \
    samtools view -N /dev/stdin $O.sam  -b > $O.bam
  samtools index $O.bam
  rm $O.sam
fi

# count aligned reads; compute cvg; get coverage stats; get split alignments

if [ ! -s $O.cvg ]    ; then cat $O.bam | bedtools bamtobed -cigar  | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai -d | tee $O.cvg  | cut -f3 | st.pl -r $S  > $O.cvg.stat ; fi
if [ ! -f $O.sa.bed ] ; then samtools view -h $O.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed ; fi

if [ ! -s $OS.00.vcf ] ; then
  cat $HP_SDIR/$HP_M.vcf > $OS.00.vcf
  echo "##sample=$S" >> $OS.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OS.00.vcf

  bcftools mpileup -f $HP_RDIR/$HP_MT.fa $O.bam -d $HP_DP | bcftools call --ploidy 2 -mv -Ov | bcftools norm -m-any  -f $HP_RDIR/$HP_MT.fa  | tee $OS.vcf | \
    fix${HP_M}Vcf.pl |  tee  $OS.fix.vcf | filterVcf.pl -sample $S -source $HP_M  |  grep -v "^#" >> $OS.00.vcf
  cat $OS.fix.vcf | maxVcf.pl | bedtools sort -header | tee $OS.max.vcf | bgzip -f -c > $OS.max.vcf.gz ; tabix -f $OS.max.vcf.gz
  annotateVcf.sh $OS.00.vcf
fi

if [ ! -s $OS.haplogroup ] ; then
  if [ "$HP_O" == "Human" ] ; then
    if [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "chrM" ] ; then  java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup
    elif [ "$HP_MT" == "RSRS" ] ; then                           java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup --rsrs
    fi
  fi
fi

if  [ ! -s $OS.fa ]  ; then
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OS.max.vcf.gz -H A | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OS.fa
  samtools faidx $OS.fa
fi

#################################################

if [ ! -s $OS.bam ] ; then
   cat $O.fa |  minimap2 $OS.fa /dev/stdin -R $RG -a | samtools view -b | samtools sort | samtools view -h > $OS.sam
   bedtools bamtobed -i $OS.sam | cut -f 1-4 | bed2bed.pl  | count.pl -i 3 -j 4  | sort | join $O.len -  -a 1 --nocheck-order | \
     perl -ane 'print "$F[0]\t$F[1]\t$F[2]\t",$F[2]/$F[1],"\n"' | perl -ane 'print  if($F[-1]>$ENV{PC});'  | cut -f1 | \
     samtools view -N /dev/stdin $OS.sam  -b > $OS.bam
  samtools index $OS.bam
  rm $OS.sam $O.fa $O.len
  touch $OS.merge.bed
fi

if [ ! -s $OSS.00.vcf ] ; then
  cat $HP_SDIR/$HP_M.vcf > $OSS.00.vcf
  echo "##sample=$S" >> $OSS.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OSS.00.vcf

  bcftools mpileup -f $OS.fa $OS.bam -d $HP_DP | bcftools call --ploidy 2 -mv -Ov | bcftools norm -m-any  -f $OS.fa  > $OSS.vcf
  cat  $OSS.vcf | \
    fix${HP_M}Vcf.pl  | \
    fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OS.max.vcf  | \
    cat $OS.max.vcf - | \
    filterVcf.pl -sample $S -source $HP_M |  bedtools sort  >> $OSS.00.vcf

  annotateVcf.sh $OSS.00.vcf
fi
