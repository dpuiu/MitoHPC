#!/usr/bin/env bash
#set -e

##############################################################################################################

# Program that setups the environmnet for the Human hs38DH reference

# Variable HP_SDIR must be pre-set !!!

##############################################################################################################

#DIRECTORY PATHS

export HP_HDIR=`readlink -f $HP_SDIR/..`	#HP home directory
export HP_BDIR=$HP_HDIR/bin/			#bin directory
export HP_JDIR=$HP_HDIR/java/			#java directory
export HP_RDIR=$HP_HDIR/RefSeq/			#reference directory

###############################################################
#SOFTWARE PATH

export PATH=$HP_SDIR:$HP_BDIR:$PATH

################################################################
#ALIGNMNET REFERENCE

#hs38DH/GRCH38(default)
export HP_RNAME=hs38DH
export HP_RMT=chrM
export HP_RNUMT="chr1:629084-634672 chr17:22521208-22521639"
export HP_RCOUNT=3366
export HP_RURL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

################################################################
#GENOME REFERENCES

export HP_O=Human		 # organism
export HP_MT=chrM                # chrM, rCRS or RSRS, FASTA file available under $HP_RDIR
export HP_MTLEN=16569
export HP_NUMT=NUMT              # NUMT FASTA file under $HP_RDIR

################################################################
#OTHER

export HP_CN=1			 # do compute copy number
export HP_L=222000               # number of MT reads to subsample; empty: no subsamling; 222000 150bp reads => ~2000x MT coverage
export HP_E=300	                 # extension(circularization)
export HP_FOPT="-q 20 -e 20"     # FASTP options: Ex: " -q 20 -e 30 "; -q: min base quality; -e: avg quality thold
export HP_DOPT="--removeDups"    # samblaster option; leave empty if no deduplication should be done
export HP_GOPT=                  # gatk mutect2 additional options : Ex "-max-reads-per-alignment-start 50" , "-mitochondria-mode"

export HP_M=mutect2 	         # SNV caller: mutect2, mutserve or freebayes
export HP_I=2		         # number of SNV iterations : 0,1,2
				 #  0: compute read counts,mtDNA-CN
                                 #  1:1 iteration (mutect2,mutserve)
                                 #  2:2 iterations (mutect2)
export HP_T1=03                  # heteroplasmy tholds
export HP_T2=05
export HP_T3=10

export HP_V=gridss               # SV caller: gridss
export HP_DP=50                 # minimum coverage

export HP_FNAME=                 # FILTERING options: Ex: noMultiallelic
export HP_FRULE=                 # FILTERING options: Ex: "grep -v multiallelic"


export HP_P=1						    # number of processors
export HP_JOPT="-Xms2G -Xmx2G -XX:ParallelGCThreads=$HP_P"  # JAVA options
export HP_MM="3G"					    # maximum memory; 4G before
################################################################
#INPUT/OUTPUT

PWD=`pwd -P`
export HP_FDIR=$PWD/fastq/      # fastq input file directory ; .fq or .fq.gz file extension
export HP_ADIR=$PWD/bams/	# bams or crams input file directory
export HP_ODIR=$PWD/out/        # output dir
export HP_IN=$PWD/in.txt        # input file to be generated

if [ -d $HP_ADIR ] ; then
  if [ ! -s $HP_IN ] ; then
    find $HP_ADIR/ -name "*.bam" -o -name "*.cram" | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
  fi
fi

###############################################################
#JOB SCHEDULING

export HP_SH="bash" ;                                                            export HP_SHS="$HP_SH"                     # bash
#export HP_SH="sbatch -J HP_$$ --cpus-per-task=$HP_P --nodes=1 --mem=5G" ;       export HP_SHS="$HP_SH -d singleton"        # SLURM
#export HP_SH="qsub -V -N HP_$$ -l mem_free=5G,h_vmem=5G -pe local $HP_P -cwd" ; export HP_SHS="$HP_SH -hold_jid HP_$$"     # SGE
