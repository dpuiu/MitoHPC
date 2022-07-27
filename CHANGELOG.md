# 2022/07/27 #

* scripts/install_sysprerequisites.sh : replaced "set -e" with "set -x"
* scripts/install_prerequisites.sh    :	replaced "set -e" with "set -x"
* RefSeq/HP.bed.gz                    : update one homopolymer interval 

# 2022/07/12 #

* scripts/install_prerequisites.sh    : updated "bwa mem" download path

# 2022/06/23 #

* scripts/filter.sh                   : commented $HP_RCOUNT check; replaced 2 "join" with "join --nocheck-order"
* scripts/intersect.pl                : added a "-header" option
* scripts/reAnnotateVcf.sh            : added the new YEALE mito scores
* scripts/annotateVcf.sh     	      :	added the new YEALE mito scores
* scripts/init.sh                     : updated HP_FRULE
* scripts/snpCount.sh                 : update "suspicious" samples 

# 2022/06/03 #

* scripts/fixmutect2secondVcf.pl      : added
* scripts/difference.pl               : added
* scripts/snpCount.sh  	       	      : added "suspicious" samples

# ... #

# 2020/11/20 #

* Added project to GitHub
