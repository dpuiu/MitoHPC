#!/usr/bin/env bash
set -ex

#######################################################################################################################################

# Program that runs the heteroplasmy pipeline on a single sample

# Input arguments
#  1: sample names
#  2: BAM/CRAM alignment file; full path
#  3: output prefix; full path

########################################################################################################################################

#set variables
export S=$1
N=`basename $2 .bam`
export N=`basename $N .cram`
IDIR=`dirname $2`
I=$IDIR/$N
ODIR=`dirname $3`; mkdir -p $ODIR
O=$3
#OA="$O.all"
ON="$O.$HP_NUMT"
OM=$O.$HP_M
OMM=$OM.$HP_M

########################################################################################################################################

#test if count and VCF output files exist; exit if they do
if [ $HP_I -eq 1 ] && [ -s $I.count ] && [ -s $OM.00.vcf  ] ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $I.count ] && [ -s $OMM.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################

#test alignment file exists and is sorted by coordinates
test -s $2
samtools view -H $2 | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

#test references match
RCOUNT=`samtools view -H $2 | grep -c "^@SQ"`
if [ $HP_RCOUNT != $RCOUNT ] ; then echo "ERROR: HP_RCOUNT=$HP_RCOUNT does not match the number of \@SQ lines=$RCOUNT in $2"; exit 1 ; fi

#generate index and indexstats files
if [ ! -s $2.bai ] && [ ! -s $2.crai ]; then samtools index $2 -@ $HP_P ; fi
if [ ! -s $I.idxstats ] ; then               samtools idxstats $2 > $I.idxstats ; fi

#get read counts
if [ ! -s $I.count ]; then cat $I.idxstats | idxstats2count.pl -sample $S -chrM $HP_RMT > $I.count ; fi

#test if VCF output files exist; exit if they do
if [ $HP_I -lt 1 ] ; then exit 0 ; fi
if [ $HP_I -eq 1 ] && [ -s $OM.00.vcf ]  ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $OMM.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################

#subsample reads
if [ ! -s $O.fq ] ; then
  R=""
  if [ $HP_L ]; then
    R=`tail -1 $I.count  | perl -ane '$R=$ENV{HP_L}/$F[-1]; print $R if($R<1)'`
    if [ $R ] ; then R="-s $R"  ; fi
  fi

  samtools view $R $2 $HP_RMT $HP_RNUMT -bu -F 0x900 -T $HP_RDIR/$HP_RNAME.fa -@ $HP_P | \
    samtools sort -n -O SAM -m $HP_MM -@ $HP_P |  \
    perl -ane 'if(/^@/) {print} elsif($P[0] eq $F[0]) {print $p,$_}; @P=@F; $p=$_;' | \
    samblaster $HP_DOPT --addMateTags | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout $HP_FOPT  > $O.fq
fi

#########################################################################################################################################

#realign subsampled reads
if  [ ! -s $O.bam ] ; then
  cat $O.fq | \
    bwa mem $HP_RDIR/$HP_MT+ - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view  -F 0x90C -h | \
    circSam.pl -ref_len $HP_RDIR/$HP_MT.fa.fai | tee $O.sam  | \
    samtools view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix | \
    count.pl -i 3 -j 4  | sort > $O.score

  cat $O.fq | \
    bwa mem $HP_RDIR/$HP_NUMT - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix  | \
    count.pl -i 3 -j 4  | sort > $ON.score

  join $O.score $ON.score -a 1 | perl -ane 'next if(@F==3 and $F[2]>$F[1]);print' | \
     intersectSam.pl $O.sam - | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P  > $O.bam

  rm -f $O.sam $O.score
  samtools index $O.bam -@ $HP_P
fi

#########################################################################################################################################

#count aligned reads; compute cvg; get stats; gets split alignments
if [ ! -s $O.count ]    ; then samtools idxstats $O.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $HP_MT > $O.count ; fi
if [ ! -s $O.cvg ]      ; then cat $O.bam | bedtools bamtobed -cigar | grep "^$HP_MT" | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai -d > $O.cvg ; fi
if [ ! -s $O.cvg.stat ] ; then cat $O.cvg | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $O.cvg.stat ; fi
if [ ! -f $O.sa.bed ]   ; then samtools view -h $O.bam  -@ $HP_P | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed ; fi

#########################################################################################################################################

#compute SNP/INDELs using mutect2 or mutserve
if [ ! -s $OM.vcf ] ; then
  if [ "$HP_M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_MT.fa -I $O.bam       -O $OM.orig.vcf $HP_GOPT --native-pair-hmm-threads $HP_P
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_MT.fa -V $OM.orig.vcf -O $OM.vcf --min-reads-per-strand 2
  elif [ "$HP_M" == "mutserve" ] ; then
    if [ "$HP_MT" == "chrM" ] ||  [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "RSRS" ] ; then
      java $HP_JOPT -jar $HP_JDIR/mutserve.jar call --deletions --insertions --level 0.01 --output $OM.vcf --reference $HP_RDIR/$HP_MT.fa $O.bam
    else
      echo "Wrong mutserve reference"
      exit 1
    fi
  elif [ "$HP_M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction 0.01 $O.bam -f $HP_RDIR/$HP_MT.fa  > $OM.vcf
  else
    echo "Unsuported SNV caller"
    exit 1
  fi
fi
rm -f $O.bam* $OM.*idx $OM.*tsv $OM.*stats

if [ ! -s $OM.00.vcf ] ; then
  # filter SNPs

  bcftools norm -m - $OM.vcf | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header> $OM.fix.vcf

  cat $HP_SDIR/$HP_M.vcf > $OM.00.vcf ; echo "##sample=$S" >> $OM.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OM.00.vcf
  cat $OM.fix.vcf | filterVcf.pl -sample $S -source $HP_M | bedtools sort >> $OM.00.vcf
  cat $OM.00.vcf | maxVcf.pl | bedtools sort -header |tee $OM.max.vcf | bgzip -f -c > $OM.max.vcf.gz ; tabix -f $OM.max.vcf.gz
  annotateVcf.sh $OM.00.vcf
fi

########################################################################################################################################

# get haplogroup
if [ "$HP_O" == "Human" ] ; then

  if [ ! -s $OM.haplogroup ] ; then
    if [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "chrM" ] ; then java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OM.vcf  --format  vcf  --out $OM.haplogroup
    elif [ "$HP_MT" == "RSRS" ] ; then                          java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OM.vcf  --format  vcf  --out $OM.haplogroup --rsrs
    fi

    if [ -s $OM.haplogroup ] ; then java $HP_JOPT -jar $HP_JDIR/haplocheck.jar --out $OM.haplocheck $OM.vcf ; fi
  fi
fi

#########################################################################################################################################

#get new consensus; format reference
if  [ ! -s $OM.fa ]  ; then
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OM.max.vcf.gz | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OM.fa
  rm -f $OM.max.vcf.gz $OM.max.vcf.gz.tbi
  samtools faidx $OM.fa
  rm -f $OM.dict; java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $OM.fa --OUTPUT $OM.dict
fi

########################################################################################################################################

#realign reads; check coverage
if  [ ! -s $OM.bam ] ; then
  samtools faidx $OM.fa $S:1-$HP_E | grep -v ">" | cat $OM.fa - > $OM+.fa
  bwa index $OM+.fa -p $OM+
  cat $O.fq | \
    bwa mem $OM+ - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view  -F 0x90C -h | \
    circSam.pl -ref_len $OM.fa.fai | tee $OM.sam  | \
    samtools  view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix  | \
    count.pl -i 3 -j 4  | sort > $OM.score
  rm -f $OM+.*

  join $OM.score $ON.score -a 1 | perl -ane 'next if(@F==3 and $F[2]>$F[1]);print'  | \
     intersectSam.pl $OM.sam - | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P  > $OM.bam

  samtools index $OM.bam -@ $HP_P
  bedtools bamtobed -i $OM.bam -ed | perl -ane 'print if($F[-2]==0);' | bedtools merge -d -3 | bed2bed.pl -min 3 > $OM.merge.bed

  rm -f $OM.sam $OM.score $ON.score $O.fq
fi


#exit if the number of iterations is set to 1
if [ $HP_I -lt 2 ] || [ $HP_M == "mutserve" ] ; then
  rm -f $OM.bam* $ON.score
  exit 0
fi

#########################################################################################################################################

#count aligned reads; compute cvg; get stats; gets split alignments
if [ ! -s $OM.count ]    ; then samtools idxstats $OM.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $S > $OM.count ; fi
if [ ! -s $OM.cvg ]      ; then cat $OM.bam | bedtools bamtobed -cigar | grep "^$S" | bedtools genomecov -i - -g $OM.fa.fai -d > $OM.cvg ; fi
if [ ! -s $OM.cvg.stat ] ; then cat $OM.cvg | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $OM.cvg.stat ; fi
if [ ! -f $OM.sa.bed ]   ; then samtools view -h $OM.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $OM.sa.bed ; fi

#########################################################################################################################################

#compute SNP/INDELs using mutect2 or mutserve
if [ ! -s $OMM.vcf ] ; then
  if [ "$HP_M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $OM.fa -I $OM.bam       -O $OMM.orig.vcf  $HP_GOPT --native-pair-hmm-threads $HP_P
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $OM.fa -V $OMM.orig.vcf -O $OMM.vcf  --min-reads-per-strand 2
  elif [ "$HP_M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction 0.01 $OM.bam -f $OM.fa  > $OMM.vcf
  else
    echo "Unsuported SNV caller"
    exit 1
  fi
fi
rm -f $OM.bam* $OMM.*idx $OMM.*tsv $OMM.*stats

if [ ! -s $OMM.00.vcf ] ; then
  # filter SNPs
  cat $HP_SDIR/$HP_M.vcf > $OMM.00.vcf ; echo "##sample=$S" >> $OMM.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OMM.00.vcf
  bcftools norm -m - $OMM.vcf | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header> $OMM.fix.vcf 
  cat $OMM.fix.vcf | fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OM.max.vcf -ffile $OM.fix.vcf | filterVcf.pl -sample $S -source $HP_M | bedtools sort  >> $OMM.00.vcf
  annotateVcf.sh $OMM.00.vcf

  # new: 2022/03/09
  intersectVcf.pl $OM.00.vcf $OM.max.vcf | cat - $OMM.00.vcf | uniqVcf.pl | bedtools sort -header > $OMM.00.vcf.tmp
  mv $OMM.00.vcf.tmp $OMM.00.vcf
fi
