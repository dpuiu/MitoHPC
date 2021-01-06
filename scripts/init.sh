#!/bin/bash -e

SDIR=`dirname $0`        # script directory
export SDIR=`readlink -f $SDIR`
export HDIR=`readlink -f $SDIR/..`
export BDIR=`readlink -f $SDIR/../bin`
export JDIR=`readlink -f $SDIR/../java`
export RDIR=`readlink -f $SDIR/../RefSeq`

export PATH=$SDIR:$BDIR:$PATH

export SH="sh"                 # bash, sbatch(SLURM), qsub(SGE,PBS)

