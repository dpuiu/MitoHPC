#!/usr/bin/bash

export IN=${1:-$PWD/in.txt}           # input file with .bam/.cram file path; required
export ODIR=${2:-out}            # output dir; should be empty; required
export M=${3:-mutect2}           # or mutserve
export H=${4:-hs38DH}            # human reference
export R=${5:-rCRS}              # or RSRS
export T1=${6:-03}               # Heteroplasmy thold
export T2=${7:-05}
export T3=${8:-10}

###############################################################

#Program that generates the heteroplasmy pipleline script

##############################################################

export SDIR=`dirname $0`        # script directory
export HDIR=`readlink -f $SDIR/..`
export RDIR=$HDIR/RefSeq
export BDIR=$HDIR/bin
export JDIR=$HDIR/java
export RDIR=$HDIR/RefSeq
export LDIR=$HDIR/lib/perl5
export SH="bash"
export SHS="bash"
export I=2

##############################################################

if [ ! -s ./init.sh ] ; then cp $SDIR/init.sh . ; fi
source ./init.sh

test -s $IN
mkdir -p $ODIR

###############################################################
printf "#!/usr/bin/bash -eux\n\n"
printf "export SDIR=$SDIR\n"
printf "export BDIR=$BDIR\n"
printf "export JDIR=$JDIR\n"
printf "export RDIR=$RDIR\n"

printf "export R=$R\n"
printf "export H=$H\n"
printf "export MT=$MT\n"
printf "export NUMT=\"$NUMT\"\n"
printf "export E=$E\n"
printf "export L=$L\n"

printf "export M=$M\n"
printf "export T1=$T1\n"
printf "export T2=$T2\n"
printf "export T3=$T3\n"

printf "export FNAME=\"$FNAME\"\n"
printf "export FRULE=\"$FRULE\"\n"

printf "export QM=$QM\n"
printf "export QA=$QA\n"

printf "export PATH=$SDIR:$BDIR:\$PATH\n"
printf "export PERLLIB=$LDIR:$PERLLIB\n"
printf "export PERL5LIB=$LDIR:$PERL5LIB\n"

printf "\n######################################\n\n"
cat $IN | perl -ane 'next if(/^#/ or @F<3);  print "$ENV{SH} $ENV{SDIR}/filter.sh $F[1] $F[2] $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{RDIR}/$ENV{R}\n";'
printf "$SHS $SDIR/getSummary.sh $IN $ODIR $M $T1 $T2 $T3\n"
##################################################################

if [ "$M" == "mutect2" ] && [ "$I" == "2" ] ; then
  printf "\n######################################\n\n"
  cat $IN | perl -ane 'next if(/^#/ or @F<3); print "$ENV{SH} $ENV{SDIR}/filter.sh $F[1] $F[2].$ENV{M} $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $F[2].$ENV{M}\n";'
  printf "$SHS $SDIR/getSummary.sh $IN $ODIR $M.$M $T1 $T2 $T3\n"
fi
