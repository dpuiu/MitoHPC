#!/usr/bin/env bash
set -ex

export M=$1 # source (mutect2...)
export T=$2 # thold

MM=$M.$M
MMM=${M}_${MM}

if [ "$#" -lt 2 ]; then exit 0 ; fi

test -s $MM.$T.concat.vcf
test -s $M.$T.concat.vcf

cat $MM.$T.concat.vcf  > $MMM.$T.concat.vcf
grep ":1$" $M.$T.concat.vcf  >> $MMM.$T.concat.vcf
snpSort.sh $MMM.$T.concat
snpCount.sh $MMM $T 1
