#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that filters SNV's based on their AF

##############################################################################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi

export M=$1 # source (mutect2...)
export T=$2 # thold

#######################################################

if [ "$#" -lt 2 ]; then exit 0 ; fi
if [ -s $HP_ODIR/$M.00.concat.vcf ] ; then
  cat $HP_ODIR/$M.00.concat.vcf | filterVcf.pl -p 0.$T > $HP_ODIR/$M.$T.concat.vcf
fi

