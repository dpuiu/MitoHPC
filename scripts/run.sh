#!/usr/bin/env bash
set -e

test -s $HP_IN
mkdir -p $HP_ODIR

###############################################################
printf "#!/usr/bin/bash -eux\n\n"
printf "export HP_SDIR=$HP_SDIR\n"
printf "export HP_BDIR=$HP_BDIR\n"
printf "export HP_JDIR=$HP_JDIR\n"
printf "export HP_RDIR=$HP_RDIR\n"
printf "export HP_ODIR=$HP_ODIR\n"

printf "export HP_R=$HP_R\n"
printf "export HP_H=$HP_H\n"
printf "export HP_MT=$HP_MT\n"
printf "export HP_NUMT=\"$HP_NUMT\"\n"
printf "export HP_E=$HP_E\n"
printf "export HP_L=$HP_L\n"

printf "export HP_M=$HP_M\n"
printf "export HP_T1=$HP_T1\n"
printf "export HP_T2=$HP_T2\n"
printf "export HP_T3=$HP_T3\n"

printf "export HP_FNAME=\"$HP_FNAME\"\n"
printf "export HP_FRULE=\"$HP_FRULE\"\n"

printf "export HP_FOPT=\"$HP_FOPT\"\n"
printf "export HP_JOPT=\"$HP_JOPT\"\n"

printf "export PATH=$HP_SDIR:$HP_BDIR:\$PATH\n"
printf "export PERLLIB=$HP_LDIR:$PERLLIB\n"
printf "export PERL5LIB=$HP_LDIR:$PERL5LIB\n"

printf "\n"
#cat $HP_IN | perl -ane 'next if(/^#/ or @F<3); $ODIR=`dirname $F[2]`; chomp $ODIR ; print "mkdir -p $ODIR ; " if($ODIR);  print "$ENV{HP_SH} $ENV{HP_SDIR}/filter$ENV{HP_I}.sh $F[1] $F[2]\n";'
cat $HP_IN | perl -ane 'next if(/^#/ or @F<3); print "$ENV{HP_SH} $ENV{HP_SDIR}/filter$ENV{HP_I}.sh $F[1] $F[2]\n";'

printf "\n"
printf "$HP_SHS $HP_SDIR/getSummary$HP_I.sh\n"

