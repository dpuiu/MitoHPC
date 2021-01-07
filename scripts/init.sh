#!/bin/bash -e

###############################################################

#Program that setups the environmnet
#Variable SDIR must be set
#Run using "source init.sh"

############################################################### 

#SDIR=`dirname $0`        # script directory
export SDIR=`readlink -f $SDIR`
export HDIR=`readlink -f $SDIR/..`
export BDIR=$HDIR/bin
export JDIR=$HDIR/java
export RDIR=$HDIR/RefSeq

export PATH=$SDIR:$BDIR:$PATH
export SH="bash"                 # bash, sbatch(SLURM), qsub(SGE,PBS)

