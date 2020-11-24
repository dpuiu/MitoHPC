#!/bin/bash -eux 
#SBATCH --job-name=bam-index
#SBATCH --partition=shared
#SBATCH --time=12:0:0
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G

P=8

I=$1
export N=`basename $I .bam` 
export N=`basename $N .cram`
export D=`dirname $I`

test -f $I
if [ ! -s $I.bai ] && [ ! -s $I.crai ]; then
  samtools index -@ $P $I
fi

if [ ! -s $D/$N.count ] ; then
  samtools view -@ $P $I -F 0x900 -c  | awk '{print $1,"all"}'    >  $D/$N.count
  samtools view -@ $P $I -F 0x904 -c  | awk '{print $1,"mapped"}' >> $D/$N.count
fi 
