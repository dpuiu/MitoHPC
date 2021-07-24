#!/usr/bin/env bash
set -e

filter.sh $1 $2       $HP_M $HP_RDIR/$HP_H $HP_RDIR/$HP_R $HP_RDIR/$HP_R
filter.sh $1 $2.$HP_M $HP_M $HP_RDIR/$HP_H $HP_RDIR/$HP_R $2.$HP_M
