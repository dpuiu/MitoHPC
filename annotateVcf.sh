#!/bin/bash -eux
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH --time=1:0:0
#SBATCH --ntasks-per-node=1

Q=$1	#sample(SRR...)
O=$2	#source(mutect2,...)
M=$3	#chrM
N=${Q%?} # new

Q=$Q.$O

HDIR=$(dirname $(readlink -f $0))/..
RDIR=$HDIR/RefSeq/

test -f $Q.vcf

#if [ ! -s $Q.00.vcf ] ; then
  cat $HDIR/scripts/$O.vcf $HDIR/scripts/$M.vcf $HDIR/scripts/header.vcf > $Q.00.vcf
  cat ${Q}.vcf | bcftools norm -m - | filterVcf.pl -sample $N -source $O | grep -v ^# | sort -k2,2n -k4,4 -k5,5 | fix${O}Vcf.pl -file $RDIR/$M.fa >> $Q.00.vcf  

  vcf-validator $Q.00.vcf
  #cat $Q.00.vcf | filterVcf.pl -p 0.01  | tee $Q.01.vcf | filterVcf.pl -p 0.03  | tee $Q.03.vcf | filterVcf.pl -p 0.05  | tee $Q.05.vcf | filterVcf.pl -p 0.10  > $Q.10.vcf

  cat $Q.00.vcf | maxVcf.pl | tee $Q.max.vcf | bgzip -f -c > $Q.max.vcf.gz  ; tabix -f $Q.max.vcf.gz
  bcftools consensus -f $RDIR/$M.fa $Q.max.vcf.gz | sed  "s|^>chrM|>$N|" | sed 's|N||' > $Q.fa

  rm $Q.mutect2.max.vcf.gz*
#fi
