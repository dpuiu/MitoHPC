#!/usr/bin/env bash
set -e

##############################################################################################################

# Program that generates the pipline script to run

##############################################################################################################

#test -s init.sh
#. ./init.sh

test -s $HP_IN

#check SNV and SV callers
if [ $HP_I != "0" ]       && [ $HP_I != "1" ]        && [ $HP_I != "2" ] ;         then  echo "Incorrect HP_I parameter" > /dev/stderr ; exit 1 ; fi
if [ $HP_M != "mutect2" ] && [ $HP_M != "mutserve" ] && [ $HP_M != "freebayes" ] ; then  echo "Incorrect HP_M parameter" > /dev/stderr ; exit 1 ; fi
if [ $HP_V ] && [ $HP_V != "gridss" ] ;                                            then  echo "Incorrect HP_V parameter" > /dev/stderr ; exit 1 ; fi

#check THOLDs
perl -e '!$ENV{HP_T1} or $ENV{HP_T1}=~/^\d\d$/ or die "ERROR : HP_T1 between 00 and 99\n"'
perl -e '!$ENV{HP_T2} or $ENV{HP_T2}=~/^\d\d$/ or die "ERROR : HP_T2 between 00 and 99\n"'
perl -e '!$ENV{HP_T3} or $ENV{HP_T3}=~/^\d\d$/ or die "ERROR : HP_T3 between 00 and 99\n"'

#test -w $HP_ADIR
mkdir -p $HP_ODIR
test -w $HP_ODIR

###############################################################
printf "#!/usr/bin/env bash\n"
printf "set -eux\n\n"

printf "export HP_SDIR=$HP_SDIR\n"
printf "export HP_BDIR=$HP_BDIR\n"
printf "export HP_JDIR=$HP_JDIR\n"
printf "export HP_RDIR=$HP_RDIR\n"
printf "export HP_ODIR=$HP_ODIR\n"

printf "export HP_RNAME=$HP_RNAME\n"
printf "export HP_RMT=$HP_RMT\n"
printf "export HP_RNUMT=\"$HP_RNUMT\"\n"
printf "export HP_RCOUNT=$HP_RCOUNT\n"
printf "export HP_RURL=\"$HP_RURL\"\n"

printf "export HP_O=$HP_O\n"
printf "export HP_MT=$HP_MT\n"
printf "export HP_MTR=$HP_MT.r$HP_E\n"
printf "export HP_MTC=$HP_MT.c$HP_E\n"
printf "export HP_MTLEN=$HP_MTLEN\n"
printf "export HP_NUMT=$HP_NUMT\n"

printf "export HP_CN=$HP_CN\n"
printf "export HP_E=$HP_E\n"
printf "export HP_L=$HP_L\n"

printf "export HP_IN=$HP_IN\n"
printf "export HP_ODIR=$HP_ODIR\n"

printf "export HP_M=$HP_M\n"
printf "export HP_I=$HP_I\n"

printf "export HP_T1=$HP_T1\n"
printf "export HP_T2=$HP_T2\n"
printf "export HP_T3=$HP_T3\n"

printf "export HP_V=$HP_V\n"
printf "export HP_DP=$HP_DP\n"

printf "export HP_FNAME=\"$HP_FNAME\"\n"
printf "export HP_FRULE=\"$HP_FRULE\"\n"

printf "export HP_FOPT=\"$HP_FOPT\"\n"
printf "export HP_DOPT=\"$HP_DOPT\"\n"
printf "export HP_GOPT=\"$HP_GOPT\"\n"
printf "export HP_JOPT=\"$HP_JOPT\"\n"
printf "export HP_MM=\"$HP_MM\"\n"
printf "export HP_P=$HP_P\n"

printf "export PATH=$HP_SDIR:$HP_BDIR:\$PATH\n"

printf "\n"
grep -v "^#" $HP_IN | sed "s|^|$HP_SH $HP_SDIR/filter.sh |"

printf "\n"
printf "$HP_SHS $HP_SDIR/getSummary.sh $HP_ODIR\n"

