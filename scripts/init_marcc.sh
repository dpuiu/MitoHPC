#!/usr/bin/env bash
#set -e

#########################################################

#Program that setups the environment on MARCC cluster (only)

##########################################################

module load gcc
module load java
module load bedtools
module load bwa
module load samtools
module load bcftools

export HP_SDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/scripts/    # script directory
export HP_JDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/java/       # java jar directory
export HP_RDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/RefSeq/     # RefSeq directory ; contains hs38DH.fa, chrM.fa, rCRS.fa ...
export HP_BDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/bin/        # executables
export HP_LDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/lib/perl5/

export PATH=$SDIR:$BDIR:$PATH
export PERLLIB=$LDIR  #:$PERLLIB
export PERL5LIB=$LDIR #:$PERL5LIB
export SH="sbatch"
