#!/usr/bin/env bash
set -ex

test -s $HP_RDIR/$HP_RNAME.fa

if [ ! -s $HP_RDIR/$HP_RNAME.bwt ] ; then
  bwa index $HP_RDIR/$HP_RNAME.fa  -p $HP_RDIR/$HP_RNAME
fi
