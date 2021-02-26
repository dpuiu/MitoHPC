#!/usr/bin/env bash 

##########################################################

#Program that summarizes aligned read counts
#Input arguments

D=$1        # input directory

##########################################################

find $D/*/ -name "*.count" | xargs cat | sort -u > $D/count.tab
