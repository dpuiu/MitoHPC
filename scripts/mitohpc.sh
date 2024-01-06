#!/usr/bin/env bash
set -ex

##############################################################################################################

# Script which launches the MitoHPC pipeline

###############################################################################################################

. $HP_SDIR/init.sh
find $HP_ADIR/  -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
$HP_SDIR/run.sh | bash
