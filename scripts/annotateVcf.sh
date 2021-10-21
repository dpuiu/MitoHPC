#!/usr/bin/env bash
set -e

#########################################################

#Program that annotates a vcf file
#Input arguments

O=$1   # input/output vcf file

##########################################################

test  -s $O

bgzip -f $O ;tabix -f $O.gz
bcftools annotate -a $HP_RDIR/HV.bed.gz    $O.gz -c "CHROM,FROM,TO,HV"  -h <(echo '##INFO=<ID=HV,Number=0,Type=Flag,Description="Hypervariable">') > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/HP.bed.gz    $O.gz -c "CHROM,FROM,TO,HP"  -h <(echo '##INFO=<ID=HP,Number=0,Type=Flag,Description="Homoloplymer">')  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/HS.bed.gz    $O.gz -c "CHROM,FROM,TO,HS"  -h <(echo '##INFO=<ID=HS,Number=0,Type=Flag,Description="Hot spot">')      > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/CDS.bed.gz   $O.gz -c "CHROM,FROM,TO,CDS"  -h <(echo '##INFO=<ID=CDS,Number=0,Type=Flag,Description="CDS">')         > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/RNR.bed.gz   $O.gz -c "CHROM,FROM,TO,RNR"  -h <(echo '##INFO=<ID=CDS,Number=0,Type=Flag,Description="RNR">')         > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/TRN.bed.gz   $O.gz -c "CHROM,FROM,TO,TRN"  -h <(echo '##INFO=<ID=CDS,Number=0,Type=Flag,Description="TRN">')         > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/DLOOP.bed.gz   $O.gz -c "CHROM,FROM,TO,DLOOP"  -h <(echo '##INFO=<ID=CDS,Number=0,Type=Flag,Description="DLOOP">')   > $O ; bgzip -f $O ; tabix -f $O.gz

bcftools annotate -a $HP_RDIR/HG.vcf.gz    $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/NUMT.vcf.gz  $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/dbSNP.vcf.gz $O.gz -c "ID"    > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/CADD.vcf.gz  $O.gz -c "INFO"  > $O 

rm -f $O.gz $O.gz.tbi
