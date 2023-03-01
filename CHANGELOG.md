# 2023/03/01 #

* updated getCN.pl, idxstats2count.pl, reAnnotateVcf.sh (slight change in the output format)
* updated the Yale MLC and MCC score to the latest ones ( received from Nicole Lake on Jan 30 2023)
* fixed getSummary.sh bug (L76: replaced max.vcf.pl with max.vcf)

# 2023/02/22 #

* updated concat2cat.pl which failed if no annotation was present

# 2023/02/06 #

* updated circFasta.sh & rotateFasta.sh 

# 2022/11/02 #

* HP_ADIR will  store the bam/cram/bai/crai files and does not have to be writable.  Optionally, it can store .idxstats & .count files. 
* If HP_ADIR/.idxstats & HP_ADIR/.count files exist, they will be copied to the sample output directoies;  otherwise they will be generated under the sample output directories;
* SNV detection on the MT 5',3' 11bp ends ; (mutect2 "bug")
* Prerequisite software updates; 
* SNV label updates (HP replaced with Homopolymer; HV with Hypervariable; HS with Hotspot)
* Removed HP_FNAME; HP_FRULE update; HP_FRULE is always applied
* HP_DP (min DP) set by default to 100

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
