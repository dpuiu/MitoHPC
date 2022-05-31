#!/bin/bash -eux

# single sample
export sample=$1

echo "project:$project"
echo "cram_path:$cram_path"
echo "bam_path:$bam_path"

###################################

dx cd
dx cd "$cram_path"
dx ls $sample.cram
dx ls $sample.cram.crai

dx mkdir -p "$bam_path"

dx run swiss-army-knife \
      -iin="$sample.cram" \
      -iin="$sample.cram.crai" \
      -icmd='samtools view -b ${in_name[0]} chr1:629084-634672 chr17:22521208-22521639 chrM -o ${in_prefix[0]}.bam' \
      --destination="$project:$bam_path" \
      -y

