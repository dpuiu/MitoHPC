#!/usr/bin/env bash 

###############################################################

#Program that runs the heteroplasmy pipeline twice, 1st using rCRS(chrM) or RSRS as reference, 2nd the new consensus
#Input arguments: 

I=$1	 	     # input .bam/.cram file ; full path
O=$2                 # output prefix
M=${3:-mutect2}      # mutect2(default) or mutserve
H=${4:-hs38DH.fa}    # human reference assembly
R=${5:-rCRS.fa}      # rCRS or RSRS.fa

##################################################################

export SDIR=`dirname $0`        # script directory
source $SDIR/init.sh
#source $SDIR/init_marcc.sh

$SDIR/filter.sh $I $O    $M $H $R $R 
$SDIR/filter.sh $I $O.$M $M $H $R $O.$M.fa

