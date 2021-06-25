#!/bin/bash -eux

M=$1			# mutect2
IN=${2:-ARIC.txt}       # in.txt
T1=03
T1=05
T1=10
F=NUMT

removeFeature1.sh $M $T1 $F $IN
removeFeature1.sh $M $T2 $F $IN
removeFeature1.sh $M $T3 $F $IN
