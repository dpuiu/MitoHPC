#!/usr/bin/env bash
#set -e

###############################################################

#Program that setups the environmnet
#Variable HP_SDIR must be pre-set

###############################################################
#PATHS

export HP_HDIR=`readlink -f $HP_SDIR/..`	#HP home directory
export HP_BDIR=$HP_HDIR/bin/			#bin directory
export HP_JDIR=$HP_HDIR/java/			#java directory
export HP_RDIR=$HP_HDIR/RefSeq/			#reference file
export HP_LDIR=					#subsample directory

export PATH=$HP_SDIR:$HP_BDIR:$PATH		#PATH

################################################################
#ORIGINAL REFERENCE

export HP_H=hs38DH          # human assembly, FASTA file available under $RDIR
export HP_MT=chrM ;         export HP_NUMT="chr1:629084-634422 chr17:22521366-22521502"                #hs38DH
#export HP_MT=M;            export HP_NUMT='1:629084-634422 17:22521366-22521502'                      #grch38_1kgmaj
#export HP_MT=NC_012920.1;  export HP_NUMT='NC_000001.11:629084-634422 NC_000017.11:22521366-22521502' #GRCh38.p13
#export HP_MT=chrM;         export HP_NUMT='chr1:564465-569708 chr17:22020692-22020827'                #hg19
#export HP_MT=MT;           export HP_NUMT='1:564465-569708 17:22020692-22020827'                      #GRCh37-lite,hs37d5
################################################################
#OTHER VARIABLES

export HP_R=chrM                 # chrM, rCRS or RSRS, FASTA file available under $HP_RDIR 
export HP_N=NUMT                 # NUMT FASTA file under $HP_RDIR
export HP_L=222000               # maximum number of MT reads; 150bp reads => ~2000x MT coverage
export HP_E=300	                 # extension(circularization)

export HP_FOPT=                  # FASTP options: Ex: " -q 20 -e 30 "

export HP_M=mutect2 	         # SNV caller: mutect2 or mutserve
export HP_I=2		         # number of mutect2 iterations : 1 or 2
export HP_T1=03                  # heteroplasmy tholds
export HP_T2=05
export HP_T3=10
export HP_FNAME=                 # FILTERING options: Ex: noMultiallelic
export HP_FRULE=                 # FILTERING options: Ex: "grep -v multiallelic"

export HP_JOPT="-Xms2G -Xmx2G -XX:ParallelGCThreads=1"  # JAVA options
export HP_MM="4G"					# maximum memory
export HP_P=1						# number of processors
################################################################
#INPUT/OUTPUT

PWD=`pwd`
export HP_ADIR=$PWD/bams/	# bams or crams
export HP_ODIR=$PWD/out/        # output dir
export HP_IN=$PWD/in.txt        # input file

###############################################################
#JOB SCHEDULING

#export HP_SH="bash" ;                                                       export HP_SHS="bash"                                      # bash
export HP_SH="sbatch  -J HP_$$  --nodes=1 --mem=5G" ;                        export HP_SHS="sbatch -J HP_$$ -d singleton"              # SLURM
#export HP_SH="qsub -V -N HP_$$ -l mem_free=5G,h_vmem=5G -pe local 1 -cwd" ; export HP_SHS="qsub -V -hold_jid HP_$$ -N HP_S$$ -cwd"    # SGE
