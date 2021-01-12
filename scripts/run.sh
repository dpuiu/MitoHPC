#!/bin/bash -e

###############################################################

#Program that generates the heteroplasmy pipleline script
#Input arguments: 

export IN=$1	 	  	 # input file with .bam/.cram file path; required
export ODIR=${2:-.}	  	 # output dir; should be empty; required
export M=${3:-mutect2}           # or mutserve
export H=${4:-hs38DH.fa}         # human reference
export R=${5:-rCRS.fa}           # or RSRS.fa

##############################################################

export SDIR=`dirname $0`        # script directory
source $SDIR/init.sh
#source $SDIR/init_marcc.sh

test -s $IN
mkdir -p $ODIR ; test -d $ODIR

###############################################################

printf "#!/bin/bash -eux\n\nexport SDIR=$SDIR\nexport JDIR=$JDIR\nexport RDIR=$RDIR\nexport H=$H\nexport R=$R\nexport PATH=$SDIR:$BDIR:\$PATH\nexport PERLLIB=$LDIR:$PERLLIB\n"

printf "\n######################################\n\n"
cat $IN | perl -ane '/.+\/(.+)\./; print "$ENV{SH} $ENV{SDIR}/filter.sh  $F[0] $ENV{ODIR}/$1/$1 $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{RDIR}/$ENV{R}\n";'
printf "\n######################################\n"

printf "\nreadCount.sh $ODIR\n"
printf "cvgCount.sh $ODIR\n"
printf "snpCount.sh $ODIR $M\n"
printf "find $ODIR -name *.$M.fa | sort | xargs cat > $ODIR/$M.fa\n"
printf "join.pl $ODIR/count.tab $ODIR/$M.tab > $ODIR/count_$M.tab\n"
printf "find $ODIR -name *.$M.merge.bed | sort | xargs cat > $ODIR/$M.merge.bed\n"

##################################################################

if [ "$M" == "mutect2" ] ; then
  printf "\n######################################\n\n"
  cat $IN | perl -ane '/.+\/(.+)\./; print "$ENV{SH} $ENV{SDIR}/filter.sh  $F[0] $ENV{ODIR}/$1/$1.$ENV{M} $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{ODIR}/$1/$1.$ENV{M}.fa\n";'
  printf "\n######################################\n"

  printf "\ncp $ODIR/$M.haplogroup.tab $ODIR/$M.$M.haplogroup.tab\n";
  printf "cvgCount.sh $ODIR $M\n"
  printf "readCount.sh $ODIR $M\n"
  printf "snpCount.sh $ODIR $M.$M\n"
  printf "join.pl $ODIR/count.$M.tab $ODIR/$M.$M.tab > $ODIR/count_$M.$M.tab\n"
fi

printf "rm -f fastp.html fastp.json\n"
