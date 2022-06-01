#!/usr/bin/env bash
set -eu

##############################################################################################################

# Program that annotates a vcf file

# Input:  vcf file
# Output: vcf file

###############################################################################################################

I=$1
O=$2

cat $I |\
  annotateVcf.pl - $HP_RDIR/MMC.vcf.gz       |\
  annotateVcf.pl - $HP_RDIR/MLC.vcf.gz       |\
  bcftools annotate -a $HP_RDIR/MCC.bed.gz -c "CHROM,FROM,TO,MCC"         -h <(echo '##INFO=<ID=MCC,Number=1,Type=String,Description="missense_OEUF">') > $O

