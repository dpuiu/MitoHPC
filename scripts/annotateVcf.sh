#!/bin/bash -eux

#1: I: input file.bam path/prefix.bam
#2: O: output

I=$1  ; test -s $I
O=$2

RDIR=$SDIR/RefSeq/
M=chrM

bgzip -f -c $I > $O.gz
tabix -f $O.gz

bcftools annotate -a $RDIR/$M.HV.bed.gz    $O.gz -c "CHROM,FROM,TO,HV"  -h <(echo '##INFO=<ID=HV,Number=0,Type=Flag,Description="Hypervariable">') > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.HP.bed.gz    $O.gz -c "CHROM,FROM,TO,HP"  -h <(echo '##INFO=<ID=HP,Number=0,Type=Flag,Description="Homoloplymer">')  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.HS.bed.gz    $O.gz -c "CHROM,FROM,TO,HS"  -h <(echo '##INFO=<ID=HS,Number=0,Type=Flag,Description="Hot spot">')      > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.CDS.bed.gz   $O.gz -c "CHROM,FROM,TO,CDS"  -h <(echo '##INFO=<ID=CDS,Number=0,Type=Flag,Description="CDS">')         > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.HG.vcf.gz    $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.NUMT.vcf.gz  $O.gz -c "INFO"  > $O ; bgzip -f $O ; tabix -f $O.gz
bcftools annotate -a $RDIR/$M.dbSNP.vcf.gz $O.gz -c "ID"    > $O
rm $O.gz $O.gz.tbi
