#!/bin/bash -e

###############################################################
# might need to edit / include this section in the .bash_profile

export SDIR=`dirname $0`        # script directory
#source $SDIR/init.sh
source $SDIR/init_marcc.sh

##################################################################

export IN=$1	 	  	 # input file with .bam/.cram file path; required
export OUT=$2	  	         # output dir; should be empty; required

export H=${3:-hs38DH.fa}         # human reference
export R=${4:-rCRS.fa}           # or RSRS.fa
export M=${5:-mutect2}           # or mutserve

test -f $IN
mkdir -p $OUT ; test -d $OUT

export SH="sbatch"
###############################################################

printf "#!/bin/bash -e\n\nexport SDIR=$SDIR\nexport JDIR=$JDIR\nexport RDIR=$RDIR\nexport H=$H\nexport R=$R\nexport M=$M\n\nPATH=$SDIR:\$PATH\n\n"  
cat $IN | perl -ane '/.+\/(.+)\./;  print "$ENV{SH} $ENV{SDIR}/filter.sh $F[0] $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{OUT}/$1 $ENV{M}\n";'                          
printf "\ncd $OUT\nreadCount.sh .\nsnpCount.sh . $M\ncd -\n"                                                                                        
exit 0
 
