#!/usr/bin/env bash

###############################################################

#Program that setups the environmnet
#Variable SDIR must be set
#Run using "source init.sh"

############################################################### 

#SDIR=`dirname $0`        # script directory
export SDIR=`readlink -f $SDIR`
export HDIR=`readlink -f $SDIR/..`
export BDIR=$HDIR/bin
export JDIR=$HDIR/java
export RDIR=$HDIR/RefSeq
export LDIR=$HDIR/lib/perl5

export PATH=$SDIR:$BDIR:$PATH
export PERLLIB=$LDIR/:$PERLLIB
export PERL5LIB=$LDIR/:$PERL5LIB

export SH="bash"                                               # parallel : filter.pl jobs
export SHS="bash"		                               # singleton :  getSummary.pl job

#export SH="sbatch -J HP$$ --partition=shared -D "                 # SLURM (MARCC)
#export SHS="sbatch -J HP$$ -d singleton --partition=shared -D $ODIR "   # SLURM

#define MYSCRATCH unless defined 
#export SH="qsub  -V -N HP$$ -wd "          # SGE (JHPCE)
#export SHS="qsub -V -hold_jid HP$$ -N HPS$$ -wd $ODIR "

################################################################

export MT=chrM		#
export NUMT='chr1:629084-634422 chr17:22521366-22521502 '   # chrM + 2 selected NUMT
export L=222000    	# maximum number of MT reads; 150bp reads => ~2000x MT coverage
export E=300       	# extension(circularization)

################################################################
#IF THE FOLLOWING LINES ARE UNCOMMENTED, THESE FOLLOWING VARIABLES WILL BE USES INSTEAD OF THE DEFAULTS OR COMMAND LINE ARGUMENTS

export H=hs38DH         # human assembly, FASTA file available under $RDIR 
export R=chrM           # chrM, rCRS or RSRS, FASTA file available under $RDIR 
export T1=03		# heteroplasmy tholds
export T2=05
export T3=10
export M=mutect2	# SNV calledr: mutect2 or mutserve
export I=2		# number of mutect2 iterations : 1 or 2

#export FNAME="noNUMT"	      # SNV filter (optional)
#export FRULE="grep -v NUMT"
#export FNAME="PASS_clustered_events"
#export FRULE='grep -P "PASS\t|clustered_events\t"' 

export QM=0		# fastp quality trimming -q param [20]; use 0 if no quality provided
export QA=0		# 			 -e       [30]
