#!/bin/bash -e
export IN=$1                     # input file with .bam/.cram file path; required
export ODIR=${2:-out}            # output dir; should be empty; required

###############################################################

#Program that generates the heteroplasmy pipleline script

##############################################################

export SDIR=`dirname $0`        # script directory
test -s ./init.sh
source ./init.sh

test -s $IN
mkdir -p $ODIR

###############################################################

printf "#!/bin/bash -eux\n\n"
printf "export SDIR=$SDIR\n"
printf "export BDIR=$BDIR\n"
printf "export JDIR=$JDIR\n"
printf "export RDIR=$RDIR\n"
printf "export R=$R\n"
printf "export HG=$HG\n"
printf "export MT=$MT\n"
printf "export NUMT=\"$NUMT\"\n"
printf "export E=$E\n"

printf "export M=$M\n"
printf "export T1=$T1\n"
printf "export T2=$T2\n"
printf "export T3=$T3\n"

printf "export PATH=$SDIR:$BDIR:\$PATH\n"
printf "export PERLLIB=$LDIR:$PERLLIB\n"

printf "\n######################################\n\n"
cat $IN | perl -ane '/.+\/(.+)\./; print "$ENV{SH} $ENV{SDIR}/filter.sh  $F[0] $ENV{ODIR}/$1/$1 $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{RDIR}/$ENV{R}\n";'
printf "\n######################################\n"

printf "\nreadCount.sh $ODIR\n"
printf "cvgCount.sh $ODIR\n"
printf "snpCount.sh $ODIR $M $T1 $T2 $T3\n"
printf "snpMerge.sh $ODIR $M $T1 $T2 $T3\n"

printf "find $ODIR -name *.$M.fa | sort | xargs cat > $ODIR/$M.fa\n"
printf "join.pl $ODIR/count.tab $ODIR/$M.tab > $ODIR/count_$M.tab\n"
printf "find $ODIR -name *.$M.merge.bed | sort | xargs cat > $ODIR/$M.merge.bed\n"

##################################################################

if [ "$M" == "mutect2" ] && [ "$I" == "2" ] ; then
  printf "\n######################################\n\n"
  cat $IN | perl -ane '/.+\/(.+)\./; print "$ENV{SH} $ENV{SDIR}/filter.sh  $F[0] $ENV{ODIR}/$1/$1.$ENV{M} $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{ODIR}/$1/$1.$ENV{M}.fa\n";'
  printf "\n######################################\n"

  printf "\ncp $ODIR/$M.haplogroup.tab $ODIR/$M.$M.haplogroup.tab\n";
  printf "cvgCount.sh $ODIR $M\n"
  printf "readCount.sh $ODIR $M\n"
  printf "snpCount.sh $ODIR $M.$M $T1 $T2 $T3\n"
  printf "join.pl $ODIR/count.$M.tab $ODIR/$M.$M.tab > $ODIR/count_$M.$M.tab\n"
fi

printf "rm -f fastp.html fastp.json\n"
