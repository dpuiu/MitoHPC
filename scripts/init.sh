#!/usr/bin/env bash
#set -e

###############################################################

#Program that setups the environmnet
#Variable HP_SDIR must be pre-set

###############################################################
#PATHS

export HP_HDIR=`readlink -f $HP_SDIR/..`
export HP_BDIR=$HP_HDIR/bin/
export HP_JDIR=$HP_HDIR/java/
export HP_RDIR=$HP_HDIR/RefSeq/
export HP_LDIR=$HP_HDIR/lib/perl5/
export PATH=$HP_SDIR:$HP_BDIR:$PATH
#export PERLLIB=$HP_LDIR:$PERLLIB
#export PERL5LIB=$HP_LDIR:$PERL5LIB

################################################################
#ORIGINAL REFERENCE

export HP_MT=chrM ;         export HP_NUMT="chr1:629084-634422 chr17:22521366-22521502"                #hs38DH
#export HP_MT=M;            export HP_NUMT='1:629084-634422 17:22521366-22521502'                      #grch38_1kgmaj
#export HP_MT=NC_012920.1;  export HP_NUMT='NC_000001.11:629084-634422 NC_000017.11:22521366-22521502' #GRCh38.p13
#export HP_MT=chrM;         export HP_NUMT='chr1:564465-569708 chr17:22020692-22020827'                #hg19
#export HP_MT=MT;           export HP_NUMT='1:564465-569708 17:22020692-22020827'                      #GRCh37-lite,hs37d5
################################################################
#OTHER VARIABLES

export HP_H=hs38DH              # human assembly, FASTA file available under $RDIR
export HP_R=chrM                # chrM, rCRS or RSRS, FASTA file available under $RDIR
export HP_L=222000              # maximum number of MT reads; 150bp reads => ~2000x MT coverage
export HP_E=300	                # extension(circularization)

export HP_FOPT=                 # FASTP options: Ex: " -q 20 -e 30 "

export HP_M=mutect2 	        # SNV caller: mutect2 or mutserve
export HP_I=2		        # number of mutect2 iterations : 1 or 2
if [[ $HP_M != "mutect2" ]] ; then export HP_I=1 ; fi
export HP_T1=03                 # heteroplasmy tholds
export HP_T2=05
export HP_T3=10
export HP_FNAME=                # FILTERING options: Ex: noMultiallelic
export HP_FRULE=                # FILTERING options: Ex: "grep -v multiallelic"

export HP_JOPT=                 # JAVA options: Ex: "--Xms2g -Xmx2g"

################################################################
#INPUT/OUTPUT

PWD=`pwd`
export HP_ADIR=$PWD/bams/	# bams or crams
export HP_ODIR=$PWD/out/        # output dir
export HP_IN=in.txt             # input file

###############################################################
#JOB SCHEDULING

export HP_SH="bash" ;                          export HP_SHS="bash"                                      # bash
#export HP_SH="sbatch --export=ALL -J HP_$$" ; export HP_SHS="sbatch --export=ALL -J HP_$$ -d singleton" # SLURM ???
#export HP_SH="sbatch -J HP_$$" ;              export HP_SHS="sbatch -J HP_$$ -d singleton"              # SLURM    
#export HP_SH="qsub  -V -N HP_$$" ;            export HP_SHS="qsub -V -hold_jid HP_$$ -N HP_S$$ "        # SGE
