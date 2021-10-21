#!/usr/bin/env bash
set -ex

getSummary.sh $HP_IN $HP_ODIR $HP_M       $HP_T1 $HP_T2 $HP_T3
getSummary.sh $HP_IN $HP_ODIR $HP_M.$HP_M $HP_T1 $HP_T2 $HP_T3
