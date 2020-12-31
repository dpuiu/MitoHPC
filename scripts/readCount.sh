#!/bin/bash -eux

D=$1        # out dir

find $D/*/ -name "*.count" | xargs grep all -H    | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' | sort -u > $D/count.all.tab
find $D/*/ -name "*.count" | xargs grep mapped -H | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' | sort -u > $D/count.mapped.tab
find $D/*/ -name "*.count" | xargs grep chrM -H   | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' | sort -u > $D/count.chrM.tab
find $D/*/ -name "*.count" | xargs grep filter -H | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' | sort -u  > $D/count.filter.tab

if [ $# -eq 1 ] ; then join.pl $D/count.all.tab $D/count.mapped.tab | join.pl - $D/count.chrM.tab | join.pl - $D/count.filter.tab | sort | egrep -v "mutect2|mutserve" | perl -ane 'BEGIN { print "Run\tall\tmapped\tchrM\tfilter\n" } print;' | getCN.pl > $D/count.tab
                  else join.pl $D/count.all.tab $D/count.mapped.tab | join.pl - $D/count.chrM.tab | join.pl - $D/count.filter.tab | sort | egrep $2 | sed "s|.$2||"    | perl -ane 'BEGIN { print "Run\tall\tmapped\tchrM\tfilter\n" } print;' | getCN.pl > $D/count.$2.tab
fi

rm $D/count.{all,mapped,chrM,filter}.tab

