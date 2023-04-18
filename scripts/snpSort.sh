#!/usr/bin/env bash
set -ex
##############################################################################################################

# Program that sorts SNVs

##############################################################################################################

if [ "$#" -lt 1 ]; then exit 0 ; fi

test -s $1.vcf
cat $1.vcf | grep "^#"  > $1.srt.vcf 

#cat $1.vcf | grep -v "^#" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> $1.srt.vcf
cat $1.vcf | grep -v "^#" | LC_ALL=C sort -k1,1 -k2,2n -k4,4 -k5,5 --parallel=$HP_P -S $HP_MM >> $1.srt.vcf
mv  $1.srt.vcf $1.vcf
