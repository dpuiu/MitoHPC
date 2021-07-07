#!/usr/bin/env bash
set -e


test -s $1.vcf	 # true heteroplasmies: Ex: examples?/RefSeq/chrM.30.h.vcf
test -s $2.vcf	 # query file
#$3: Ex: 0.03     # heteroplasmiy tholds (min)
#$4: Ex: 0.97     #                      (max)

intersectVcf.pl  $2.vcf $1.vcf -sm -min2 $3 -max2 $4 > $2.tp.vcf
differenceVcf.pl $2.vcf $1.vcf -sm -min2 $3 -max2 $4 > $2.fp.vcf

