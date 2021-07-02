#!/bin/bash -eux

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
bcftools annotate -a $HP_RDIR/HG.vcf.gz    $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/NUMT.vcf.gz  $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $HP_RDIR/dbSNP.vcf.gz $O.gz -c "ID"    > $O 
rm -f $O.gz $O.gz.tbi
