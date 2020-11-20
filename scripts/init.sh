#!/bin/bash -eux

export SDIR=`dirname $0`        # script directory
export JDIR=$SDIR/../java/      # java jar directory
export RDIR=$SDIR/../RefSeq/    # RefSeq directory ; contains hs38DH.fa, chrM.fa, rCRS.fa ...

export H=hs38DH.fa              # human reference
export R=rCRS.fa                # or RSRS.fa
export M=mutect2                # or mutserve
export PATH=$SDIR:$PATH


