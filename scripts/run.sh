#!/bin/bash -e

###############################################################
# might need to edit / include this section in the .bash_profile

export SDIR=`dirname $0`        # script directory
source $SDIR/init.sh
#source $SDIR/init_marcc.sh

##################################################################

export IN=$1	 	  	 # input file with .bam/.cram file path; required
export ODIR=${2:-.}	  	 # output dir; should be empty; required
export M=${3:-mutect2}           # or mutserve
export H=${4:-hs38DH.fa}         # human reference
export R=${5:-rCRS.fa}           # or RSRS.fa

test -f $IN
mkdir -p $ODIR ; test -d $ODIR
###############################################################

printf "#!/bin/bash -eux\n\nexport SDIR=$SDIR\nexport JDIR=$JDIR\nexport RDIR=$RDIR\nexport H=$H\nexport R=$R\n\nexport PATH=$SDIR:\$PATH\n\n"

cat $IN | perl -ane '/.+\/(.+)\./; print "$ENV{SH} $ENV{SDIR}/filter.sh  $F[0] $ENV{ODIR}/$1/$1 $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{RDIR}/$ENV{R}\n";'
printf "\n"

printf "\n"
printf "readCount.sh $ODIR\n"
printf "snpCount.sh $ODIR $M\n"
printf "find $ODIR -name *.$M.fa | xargs cat > $M.fa\n"
printf "join.pl count.tab $M.tab > count_$M.tab\n"

exit 0
