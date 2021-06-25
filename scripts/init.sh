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

################################################################

export MT=chrM						    # hs38DH
export NUMT='chr1:629084-634422 chr17:22521366-22521502 '   # hs38DH

#export MT=MT						    # hg19
#export NUMT='1:564465-569708 17:22020692-22020827 '	    # hg19

################################################################
#IF THE FOLLOWING LINES ARE UNCOMMENTED, THESE FOLLOWING VARIABLES WILL BE USES INSTEAD OF THE DEFAULTS OR COMMAND LINE ARGUMENTS

export H=hs38DH         # human assembly, FASTA file available under $RDIR
export R=chrM           # chrM, rCRS or RSRS, FASTA file available under $RDIR
export L=222000         # maximum number of MT reads; 150bp reads => ~2000x MT coverage
export E=300	        # extension(circularization)
export T1=03		# heteroplasmy tholds
export T2=05
export T3=10
export M=mutect2	# SNV caller: mutect2 or mutserve
export I=2		# number of mutect2 iterations : 1 or 2

#export FNAME="noNUMT"	      # SNV filter (optional)
#export FRULE="grep -v NUMT"
#export FNAME="PASS_clustered_events"
#export FRULE='grep -P "PASS\t|clustered_events\t"' 

export FOPT=""          # FASTP options : Ex: " -q 20 -e 30 "

################################################################
#INPUT/OUTPUT

PWD=`pwd`
export ADIR="${ADIR:-$PWD/bams/}"	# bams or crams
export ODIR="${ODIR:-$PWD/out/}"  ; mkdir -p $ODIR
if [ -z "$IN" ] ; then export IN=$PWD/in.txt  ; find $ADIR/ | egrep "\.bam$|\.cram$" | $SDIR/ls2in.pl -out $ODIR | sort > $IN ; fi

###############################################################
#JOB SCHEDULING 

export SH="bash"                                                      # bash 
export SHS="bash"

#export SH="sbatch -p shared -J HP$$ -D "                             # SLURM
#export SHS="sbatch -p shared -J HP$$ -d singleton -D $ODIR "

#export SH="qsub  -V -N HP$$ -wd "                                    # SGE
#export SHS="qsub -V -hold_jid HP$$ -N HPS$$ -wd $ODIR "

