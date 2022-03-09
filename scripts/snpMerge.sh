#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that merges 1st & 2nd iteration SNV's

##############################################################################################################

if [ "$#" -lt 3 ]; then exit 0 ; fi

export M=$1   # source (mutect2...)
export MM=$2  # source (mutect2.mutect2)
export T=$3   # thold

#######################################################

test -f $HP_ODIR/$M.00.concat.vcf
test -f $HP_ODIR/$MM.00.concat.vcf

#cat $HP_ODIR/$M.00.concat.vcf  | filterVcf.pl -p 0.$T > $HP_ODIR/$M.$T.concat.vcf
#cat $HP_ODIR/$MM.00.concat.vcf | filterVcf.pl -p 0.$T > $HP_ODIR/$MM.$T.concat.vcf


cat $HP_ODIR/$M.$T.concat.vcf  | grep ":1$" | cat $HP_ODIR/$MM.$T.concat.vcf - | bedtools sort -header | uniqVcf.pl > $HP_ODIR/$MM.$T.concat.vcf.tmp
mv $HP_ODIR/$MM.$T.concat.vcf.tmp $HP_ODIR/$MM.$T.concat.vcf
