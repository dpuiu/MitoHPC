#!/bin/bash -eux

D=$1        # out dir

find $D/*/ -name "*.count" | xargs grep all    | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' > $D/count.all.tab
find $D/*/ -name "*.count" | xargs grep mapped | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' > $D/count.mapped.tab
find $D/*/ -name "*.count" | xargs grep chrM   | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' > $D/count.chrM.tab
find $D/*/ -name "*.count" | xargs grep filter | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);' > $D/count.filter.tab

join.pl $D/count.all.tab $D/count.mapped.tab | join.pl - $D/count.chrM.tab | join.pl - $D/count.filter.tab  | perl -ane 'BEGIN { print "Run\tall\tmapped\tchrM\tfilter\n" } print;' |  sort > $D/count.tab
rm $D/count.*.tab
