#!/bin/bash -eux

export SDIR=`dirname $0`        # script directory
export JDIR=$SDIR/../java/      # java jar directory
export RDIR=$SDIR/../RefSeq/    # RefSeq directory ; contains hs38DH.fa, chrM.fa, rCRS.fa ...
export PATH=$SDIR:$PATH
export SH="sh"                 # bash, sbatch(SLURM), qsub(SGE,PBS)

