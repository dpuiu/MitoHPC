#!/bin/bash -e

SDIR=`dirname $0`        # script directory

export SDIR=`readlink -f $SDIR`
export PATH=$SDIR:$SDIR/../bin/:$PATH
echo $PATH

export JDIR=`readlink -f $SDIR/../java`
export RDIR=`readlink -f $SDIR/../RefSeq`

export SH="sh"                 # bash, sbatch(SLURM), qsub(SGE,PBS)

