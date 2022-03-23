#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that generates SNV counts

##############################################################################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi

export M=$1   # source (mutect2...)
export T=$2   # thold

#######################################################

test -f $HP_ODIR/$M.$T.concat.vcf
cat $HP_ODIR/$M.$T.concat.vcf | eval $HP_FRULE > $HP_ODIR/$M.$HP_FNAME.$T.concat.vcf
