#!/bin/bash -eux

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

#export H=hs38DH.fa              # human reference
#export R=rCRS.fa                # or RSRS.fa
#export M=mutect2                # or mutserve
export PATH=$SDIR:$BDIR:$PATH

echo Success!
