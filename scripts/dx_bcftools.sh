#!/bin/bash -eux

# single sample
export sample=$1
export dir=`dirname $sample`

echo "project:$project"
echo "vcfgz_path:$vcfgz_path"
echo "vcf_path:$vcf_path"

####################################################

dx cd
dx mkdir -p "$vcf_path/$dir"

dx cd "$vcfgz_path"
dx ls $sample.g.vcf.gz
dx ls $sample.g.vcf.gz.tbi

dx run swiss-army-knife \
      -iin=$sample.g.vcf.gz \
      -iin=$sample.g.vcf.gz.tbi \
      -icmd='bcftools filter ${in_name[0]} -r chrM > ${in_prefix[0]}.chrM.vcf' \
      --destination=$project:$vcf_path/$dir \
      -y
