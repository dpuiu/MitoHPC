#!/bin/bash -eux

# single sample
export sample=$1

dx cd
dx cd "$bam_path"
dx ls "$sample.bam"

dx run swiss-army-knife \
      -iin=$sample.bam \
      -icmd='samtools index $in_name' \
      -y
