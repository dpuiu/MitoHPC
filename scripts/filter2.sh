#!/bin/bash -eux

###############################################################
# might need to edit / include this section in the .bash_profile

export SDIR=`dirname $0`        # script directory
source $SDIR/init.sh
#source $SDIR/init_marcc.sh

##################################################################

I=$1	 	  	   # input file with .bam/.cram file path; required
O=$2
M=${3:-mutect2}           # or mutserve
H=${4:-hs38DH.fa}    # human reference
R=${5:-rCRS.fa}      # or RSRS.fa

$SDIR/filter.sh $I $O    $M $H $R $R 
$SDIR/filter.sh $I $O.$M $M $H $R $O.$M.fa

