#!/bin/bash -eux

D=$1        # out dir

find $D/ -name "*.count" | xargs grep all | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);'     > count.all.tab
find $D/ -name "*.count" | xargs grep mapped | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);'  > count.mapped.tab
find $D/ -name "*.count" | xargs grep chrM | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);'    > count.chrM.tab
find $D/ -name "*.count" | xargs grep filter | perl -ane 'print "$1\t$2\n" if(/.+\/(\S+).count:(\d+)/);'  > count.filter.tab

join.pl count.all.tab count.mapped.tab | join.pl - count.chrM.tab | join.pl - count.filter.tab  | perl -ane 'BEGIN { print "Run\tall\tmapped\tchrM\tfilter\n" } print;' | column -t  | sort > count.tab
rm count.*.tab
