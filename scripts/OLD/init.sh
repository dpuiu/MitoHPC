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
export LDIR=$HDIR/lib/perl5

export PATH=$SDIR:$BDIR:$PATH
export PERLLIB=$LDIR/:$PERLLIB
export SH="bash"                 # bash, sbatch(SLURM), qsub(SGE,PBS)

################################################################

export HG=hs38DH
export MT=chrM
export R=rCRS
export NUMT='chr1:629084-634422 chr17:22521366-22521502 '   # chrM + 2 selected NUMT
export L=222000    # ~2000x MT coverage
export E=300       # extension(circularization)
export T1=03
export T2=05
export T3=10

######################################

#INPUT/OUTPUT
export IN=in.txt
export ODIR=out/
export M=mutect2
