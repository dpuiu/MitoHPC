#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that creates a bwa index

# Input:  genome FASTA file
# Output: bwa index files

##############################################################################################################

test -s $HP_RDIR/$HP_RNAME.fa

if [ ! -s $HP_RDIR/$HP_RNAME.bwt ] ; then
  bwa index $HP_RDIR/$HP_RNAME.fa  -p $HP_RDIR/$HP_RNAME
fi
