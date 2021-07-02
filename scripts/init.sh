#!/usr/bin/env bash

###############################################################

#Program that setups the environmnet
#Variable HP_SDIR must be pre-set

###############################################################

#HP_SDIR=`dirname $0`        # script directory
export HP_SDIR=`readlink -f $HP_SDIR`
export HP_HDIR=`readlink -f $HP_SDIR/..`
export HP_BDIR=$HP_HDIR/bin
export HP_JDIR=$HP_HDIR/java
export HP_RDIR=$HP_HDIR/RefSeq
export HP_LDIR=$HP_HDIR/lib/perl5

export PATH=$HP_SDIR:$HP_BDIR:$PATH
export PERLLIB=$HP_LDIR/:$PERLLIB
export PERL5LIB=$HP_LDIR/:$PERL5LIB

################################################################
#ORIGINAL REFERENCE

#hs38DH
export HP_MT="${MT:-chrM}"
export HP_NUMT="${NUMT:-chr1:629084-634422 chr17:22521366-22521502}"

#grch38_1kgmaj
#export HP_MT=M
#export HP_NUMT='1:629084-634422 17:22521366-22521502'

#GRCh38.p13
#export HP_MT=NC_012920.1
#export HP_NUMT='NC_000001.11:629084-634422 NC_000017.11:22521366-22521502'

#############
#hg19
#export HP_MT=chrM
#export HP_NUMT='chr1:564465-569708 chr17:22020692-22020827'

#GRCh37-lite,hs37d5
#export HP_MT=MT
#export HP_NUMT='1:564465-569708 17:22020692-22020827'

################################################################

export HP_H="${HP_H:-hs38DH}"              # human assembly, FASTA file available under $RDIR
export HP_R="${HP_R:-chrM}"                # chrM, rCRS or RSRS, FASTA file available under $RDIR
export HP_L="${HP_L:-222000}"              # maximum number of MT reads; 150bp reads => ~2000x MT coverage
export HP_E="${HP_E:-300}"	           # extension(circularization)
export HP_T1="${HP_T1:-03}"		   # heteroplasmy tholds
export HP_T2="${HP_T2:-05}"
export HP_T3="${HP_T3:-10}"
export HP_M="${HP_M:-mutect2}" 	           # SNV caller: mutect2 or mutserve
export HP_I="${HP_I:-2}"		   # number of mutect2 iterations : 1 or 2
export HP_FNAME="${HP_FNAME:-}"            # Ex: noNUMT ;            PASS_clustered_events
export HP_FRULE="${HP_FRULE:-}"            # Ex: 'grep -v NUMT'      'grep -P "PASS\t|clustered_events\t"'
export HP_FOPT="${HP_FOPT:-}"              # FASTP options : Ex: " -q 20 -e 30 "
export HP_JOPT="${HP_JOPT:--Xms2g -Xmx2g}"  # JAVA options

if [[ $HP_M != "mutect2" ]] ; then export HP_I=1 ; fi


################################################################
#INPUT/OUTPUT

PWD=`pwd`
export HP_ADIR="${HP_ADIR:-$PWD/bams/}"	# bams or crams
export HP_ODIR="${HP_ODIR:-$PWD/out/}"  ; mkdir -p $HP_ODIR
export HP_IN="${HP_IN:-$PWD/in.txt}"
if [ ! -s $HP_IN ] ;  then find $HP_ADIR/ -type f  -name "*.bam" -o -name "*.cram" | ls2in.pl -out $HP_ODIR | sort > $HP_IN ; fi

###############################################################
#JOB SCHEDULING

#bash
export HP_SH="${HP_SH:-bash}"
export HP_SHS="${HP_SHS:-bash}"

#SLURM
#export HP_SH="sbatch -p shared -J HP_$$ -D "
#export HP_SHS="sbatch -p shared -J HP_$$ -d singleton -D $HP_ODIR "

#SGE
#export HP_SH="qsub  -V -N HP_$$ -wd "
#export HP_SHS="qsub -V -hold_jid HP_$$ -N HP_S$$ -wd $HP_ODIR "

################################################################
#GENERATE SCRIPTS

#if [ ! -s $PWD/samtools.all.sh ] ; then
#  cut -f2 $HP_IN  | sed "s|^|$HP_SH $HP_SDIR/samtools.sh |" > samtools.all.sh
#  echo "cut -f2 $HP_IN | sed -r 's|(.*)\.|\1\t|g' | cut -f1 | sed 's|$|.count|' | xargs cat | uniq.pl | getCN.pl > $HP_ODIR/count.tab" >> samtools.all.sh
#fi

#if [ ! -s $PWD/filter.all.sh ] ; then
#   $HP_SDIR/run.sh > filter.all.sh
#fi
