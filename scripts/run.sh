#!/usr/bin/bash

test -s ./init.sh
. ./init.sh

test -s $IN
mkdir -p $ODIR

###############################################################
printf "#!/usr/bin/bash -eux\n\n"
printf "export SDIR=$SDIR\n"
printf "export BDIR=$BDIR\n"
printf "export JDIR=$JDIR\n"
printf "export RDIR=$RDIR\n"
printf "export ODIR=$ODIR\n"

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

printf "export FOPT=\"$FOPT\"\n"

printf "export PATH=$SDIR:$BDIR:\$PATH\n"
printf "export PERLLIB=$LDIR:$PERLLIB\n"
printf "export PERL5LIB=$LDIR:$PERL5LIB\n"

printf "\n######################################\n\n"
cat $IN | perl -ane 'next if(/^#/ or @F<3); $ODIR=`dirname $F[2]`; chomp $ODIR ; print "mkdir -p $ODIR ; " if($ODIR); $ODIR="" unless($ENV{SH}=~/^qsub/ or $ENV{SH}=~/^sbatch/); print "$ENV{SH} $ODIR $ENV{SDIR}/filter.sh $F[1] $F[2] $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $ENV{RDIR}/$ENV{R}\n";'
printf "$SHS $SDIR/getSummary.sh $IN $ODIR $M $T1 $T2 $T3\n"

##################################################################

if [ "$M" == "mutect2" ] && [ "$I" == "2" ] ; then
  printf "\n######################################\n\n"
  cat $IN | perl -ane 'next if(/^#/ or @F<3); $ODIR=`dirname $F[2]`; chomp $ODIR ; $ODIR="" unless($ENV{SH}=~/^qsub/ or $ENV{SH}=~/^sbatch/); print "$ENV{SH} $ODIR $ENV{SDIR}/filter.sh $F[1] $F[2].$ENV{M} $ENV{M} $ENV{RDIR}/$ENV{H} $ENV{RDIR}/$ENV{R} $F[2].$ENV{M}\n";'
  printf "$SHS $SDIR/getSummary.sh $IN $ODIR $M.$M $T1 $T2 $T3\n"
fi
