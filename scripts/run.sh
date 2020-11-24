#!/bin/bash -e

###############################################################
# might need to edit / include this section in the .bash_profile

export SDIR=`dirname $0`        # script directory
source $SDIR/init.sh
#source $SDIR/init_marcc.sh

##################################################################

export IN=$1	 	  	 # input file with .bam/.cram file path; required
export OUT=${2:-.}	  	 # output dir; should be empty; required
export H=${3:-hs38DH.fa}         # human reference
export R=${4:-rCRS.fa}           # or RSRS.fa
export M=${5:-mutect2}           # or mutserve

test -f $IN
mkdir -p $OUT ; test -d $OUT
###############################################################

printf "#!/bin/bash -eux\n\nexport SDIR=$SDIR\nexport JDIR=$JDIR\nexport RDIR=$RDIR\nexport H=$H\nexport R=$R\nexport M=$M\n\nexport PATH=$SDIR:\$PATH\n\n"
cat $IN | perl -ane '/.+\/(.+)\./;  print "$ENV{SH} $ENV{SDIR}/filter.sh $F[0] $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $1 $ENV{M}\n";'
printf "\n"
printf "readCount.sh $OUT\n"
printf "snpCount.sh $OUT $M\n"
printf "find $OUT -name *.$M.fa | xargs cat > $M.fa\n"
printf "find join.pl count.tab $M.tab > count_$M.tab\n"

exit 0
