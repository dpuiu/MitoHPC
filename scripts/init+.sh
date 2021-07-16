#!/usr/bin/env bash
set -e

mkdir -p $HP_ODIR
if [ ! -s $HP_IN ] ;
  then find $HP_ADIR/ -type f  -name "*.bam" -o -name "*.cram" | ls2in.pl -out $HP_ODIR | sort > $HP_IN
fi

################################################################
#GENERATE SCRIPTS

if [ ! -s $PWD/samtools.all.sh ] ; then
  cut -f2 $HP_IN  | sed "s|^|$HP_SH $HP_SDIR/samtools.sh |" > samtools.all.sh
  echo "cut -f2 $HP_IN | sed -r 's|(.*)\.|\1\t|g' | cut -f1 | sed 's|$|.count|' | xargs cat | uniq.pl | getCN.pl > $HP_ODIR/count.tab" >> samtools.all.sh
fi

if [ ! -s $PWD/filter.all.sh ] ; then
   $HP_SDIR/run.sh > filter.all.sh
fi
