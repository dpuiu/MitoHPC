#!/bin/bash -eux

##########################################################

#Program that summarizes aligned read counts
#Input arguments

D=$1 			# out dir
M=$2 			# snp caller: mutect2 or mutserve
T1=${3:-03}             # Heteroplasmy tholds
T2=${4:-05}
T3=${5:-10}

snpMerge1.sh $D $M $T1
snpMerge1.sh $D $M $T2
snpMerge1.sh $D $M $T3
