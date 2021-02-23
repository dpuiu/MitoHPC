#!/usr/bin/env bash 
#########################################################

#Program that setups the environment on MARCC cluster (only)

##########################################################

module load gcc
module load java
module load bedtools
module load bwa
module load samtools
module load bcftools

export SDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/scripts/   # script directory
export JDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/java/      # java jar directory
export RDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/RefSeq/    # RefSeq directory ; contains hs38DH.fa, chrM.fa, rCRS.fa ...
export BDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/bin/       # executables
export LDIR=/work-zfs/darking1/active/projects/mito/Heteroplasmy/HP/lib/perl5

export PATH=$SDIR:$BDIR:$PATH
export PERLLIB=$LDIR:$PERLLIB
export SH="sbatch"
