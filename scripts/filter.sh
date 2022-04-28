#!/usr/bin/env bash
set -e

#########################################################################################################################################
# Program that runs the heteroplasmy pipeline on a single sample

# Input arguments
#  1: sample names
#  2: BAM/CRAM alignment file; full path
#  3: output prefix; full path
#########################################################################################################################################

#set variables
export S=$1		# sample name
N=`basename $2 .bam`
export N=`basename $N .cram`
IDIR=`dirname $2`
I=$IDIR/$N
ODIR=`dirname $3`; mkdir -p $ODIR
O=$3
ON=$O.$HP_NUMT
OS=$O.$HP_M
OSS=$OS.$HP_M

#########################################################################################################################################
# test if count and VCF output files exist; exit if they do

if [ $HP_I -lt 1 ] && [ -s $I.count ]    ; then exit 0 ; fi
if [ $HP_I -eq 1 ] && [ -s $OS.00.vcf  ] ; then exit 0 ; fi
if [ $HP_I -ge 2 ] && [ -s $OSS.00.vcf ] ; then exit 0 ; fi

#########################################################################################################################################
# test alignment file exists and is sorted by coordinates

# test alignment input file exists and is sorted 
test -s $2
samtools view -H $2 | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

# test references match
RCOUNT=`samtools view -H $2 | grep -c "^@SQ"`
if [ $HP_RCOUNT != $RCOUNT ] ; then echo "ERROR: HP_RCOUNT=$HP_RCOUNT does not match the number of \@SQ lines=$RCOUNT in $2"; exit 1 ; fi

# generate index and indexstats files
if [ ! -s $2.bai ] && [ ! -s $2.crai ]; then samtools index $2 -@ $HP_P ; fi
if [ ! -s $I.idxstats ] ; then               samtools idxstats $2 > $I.idxstats ; fi

# get read counts
if [ ! -s $I.count ]; then cat $I.idxstats | idxstats2count.pl -sample $S -chrM $HP_RMT > $I.count ; fi

# test if there are any MT reads
MTCOUNT=`tail -1 $I.count| cut -f4`
if [ $MTCOUNT -lt 1 ]; then echo "ERROR: There are no MT reads in $2; plese remove it from $HP_IN" ; exit 1 ; fi

#########################################################################################################################################
# subsample reads

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
# realign subsampled reads

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

  # new April 2022
  #join $O.score $ON.score -a 1 | perl -ane 'print if(@F==3 and $F[2]>$F[1])' | \
  #   tee $ON.score1 | \
  #   intersectSam.pl $O.sam - | \
  #   samtools view -bu | \
  #   samtools sort -m $HP_MM -@ $HP_P  > $ON.bam

  #rm -f $O.sam $O.score
  samtools index $O.bam  -@ $HP_P
  #samtools index $ON.bam -@ $HP_P
fi
#########################################################################################################################################
# count aligned reads; compute cvg; get coverage stats; get split alignments

if [ ! -s $O.count ]    ; then samtools idxstats $O.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $HP_MT > $O.count ; fi
if [ ! -s $O.cvg ]      ; then cat $O.bam | bedtools bamtobed -cigar | grep "^$HP_MT" | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai -d > $O.cvg ; fi
if [ ! -s $O.cvg.stat ] ; then cat $O.cvg | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $O.cvg.stat ; fi
if [ ! -f $O.sa.bed ]   ; then samtools view -h $O.bam  -@ $HP_P | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed ; fi

#########################################################################################################################################
# compute SNVs using mutect2/mutserve/freebayes

if [ ! -s $OS.vcf ] ; then
  if [ "$HP_M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_MT.fa -I $O.bam       -O $OS.orig.vcf $HP_GOPT --native-pair-hmm-threads $HP_P
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_MT.fa -V $OS.orig.vcf -O $OS.vcf --min-reads-per-strand 2

    #new april 2022
    #java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_NUMT.fa -I $ON.bam      -O $ON.orig.vcf $HP_GOPT --native-pair-hmm-threads $HP_P
    #java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_NUMT.fa -V $ON.orig.vcf -O $ON.vcf --min-reads-per-strand 2

  elif [ "$HP_M" == "mutserve" ] ; then
    if [ "$HP_MT" == "chrM" ] ||  [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "RSRS" ] ; then
      java $HP_JOPT -jar $HP_JDIR/mutserve.jar call --deletions --insertions --level 0.01 --output $OS.vcf --reference $HP_RDIR/$HP_MT.fa $O.bam
    else
      echo "Wrong mutserve reference"
      exit 1
    fi
  elif [ "$HP_M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction 0.01 $O.bam -f $HP_RDIR/$HP_MT.fa  > $OS.vcf
  else
    echo "Unsuported SNV caller"
    exit 1
  fi
fi

if [ ! -s $OS.00.vcf ] ; then
  # filter SNVs
  bcftools norm -m - $OS.vcf | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header> $OS.fix.vcf
  cat $HP_SDIR/$HP_M.vcf > $OS.00.vcf ; echo "##sample=$S" >> $OS.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OS.00.vcf
  cat $OS.fix.vcf | filterVcf.pl -sample $S -source $HP_M | bedtools sort >> $OS.00.vcf  # to add -depth $HP_DP
  cat $OS.00.vcf | maxVcf.pl | bedtools sort -header |tee $OS.max.vcf | bgzip -f -c > $OS.max.vcf.gz ; tabix -f $OS.max.vcf.gz
  annotateVcf.sh $OS.00.vcf
fi

#########################################################################################################################################
#  identify SVs


if [ $HP_V ] && [ "$HP_V" == "gridss" ] ; then
  OV=$O.$HP_V
  if [ ! -s $OV.00.vcf ] ; then
    gridss --jar $HP_JDIR/gridss.jar -r $HP_RDIR/$HP_MT.fa -o $OV.vcf.gz $O.bam  -t 1 -w $OV
    cat $HP_SDIR/gridss.vcf > $OV.fix.vcf
    bcftools view -i 'FILTER="PASS"' $OV.vcf.gz | bcftools query  -f "%CHROM\t%POS\t.\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT:DP:AF\t[%GT:%REF:%AF]\n" | grep -v -P 'h\t'  >> $OV.fix.vcf
    cat $OV.fix.vcf | filterVcf.pl -sample $S -source $HP_V  > $OV.00.vcf #  to add -depth $HP_DP
    rm -rf $OV $OV.vcf.gz.* $OV.fix.vcf
  fi
fi

########################################################################################################################################
# get haplogroup

if [ "$HP_O" == "Human" ] ; then

  if [ ! -s $OS.haplogroup ] ; then
    if [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "chrM" ] ; then java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup
    elif [ "$HP_MT" == "RSRS" ] ; then                          java $HP_JOPT -jar $HP_JDIR/haplogrep.jar classify --in $OS.vcf  --format  vcf  --out $OS.haplogroup --rsrs
    fi

    if [ -s $OS.haplogroup ] ; then java $HP_JOPT -jar $HP_JDIR/haplocheck.jar --out $OS.haplocheck $OS.vcf ; fi
  fi
fi

#########################################################################################################################################
# get new consensus; format reference

if  [ ! -s $OS.fa ]  ; then
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OS.max.vcf.gz | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OS.fa
  rm -f $OS.max.vcf.gz $OS.max.vcf.gz.tbi
  samtools faidx $OS.fa
  rm -f $OS.dict; java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $OS.fa --OUTPUT $OS.dict
fi

########################################################################################################################################
# realign reads; check coverage

if  [ ! -s $OS.bam ] ; then
  samtools faidx $OS.fa $S:1-$HP_E | grep -v ">" | cat $OS.fa - > $OS+.fa
  bwa index $OS+.fa -p $OS+
  cat $O.fq | \
    bwa mem $OS+ - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" | \
    samtools view  -F 0x90C -h | \
    circSam.pl -ref_len $OS.fa.fai | tee $OS.sam  | \
    samtools  view -bu | \
    bedtools bamtobed -i /dev/stdin -tag AS | bed2bed.pl -rmsuffix  | \
    count.pl -i 3 -j 4  | sort > $OS.score
  rm -f $OS+.*

  join $OS.score $ON.score -a 1 | perl -ane 'next if(@F==3 and $F[2]>$F[1]);print'  | \
     intersectSam.pl $OS.sam - | \
     samtools view -bu | \
     samtools sort -m $HP_MM -@ $HP_P  > $OS.bam

  samtools index $OS.bam -@ $HP_P
  bedtools bamtobed -i $OS.bam -ed | perl -ane 'print if($F[-2]==0);' | bedtools merge -d -3 | bed2bed.pl -min 3 > $OS.merge.bed

  rm -f $OS.sam $OS.score $ON.score $O.fq
fi

rm -f $O.bam*
rm -f $OS.*idx $OS.*tsv $OS.*stats

# exit if the number of iterations is set to 1
if [ $HP_I -lt 2 ] || [ $HP_M == "mutserve" ] ; then
  rm -f $OS.bam* $ON.score
  exit 0
fi

#########################################################################################################################################
# 2nd iteration
# count aligned reads; compute cvg; compute cvg stats; gets split alignments

if [ ! -s $OS.count ]    ; then samtools idxstats $OS.bam  -@ $HP_P | idxstats2count.pl -sample $S -chrM $S > $OS.count ; fi
if [ ! -s $OS.cvg ]      ; then cat $OS.bam | bedtools bamtobed -cigar | grep "^$S" | bedtools genomecov -i - -g $OS.fa.fai -d > $OS.cvg ; fi
if [ ! -s $OS.cvg.stat ] ; then cat $OS.cvg | cut -f3 | st.pl | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $OS.cvg.stat ; fi
if [ ! -f $OS.sa.bed ]   ; then samtools view -h $OS.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $OS.sa.bed ; fi

#########################################################################################################################################
# compute SNP/INDELs using mutect2/freebayes
if [ ! -s $OSS.vcf ] ; then
  if [ "$HP_M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $OS.fa -I $OS.bam       -O $OSS.orig.vcf  $HP_GOPT --native-pair-hmm-threads $HP_P
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $OS.fa -V $OSS.orig.vcf -O $OSS.vcf  --min-reads-per-strand 2
  elif [ "$HP_M" == "freebayes" ] ; then
    freebayes -p 1 --pooled-continuous --min-alternate-fraction 0.01 $OS.bam -f $OS.fa  > $OSS.vcf
  else
    echo "Unsuported SNV caller"
    exit 1
  fi
fi

if [ ! -s $OSS.00.vcf ] ; then
  # filter SNPs
  cat $HP_SDIR/$HP_M.vcf > $OSS.00.vcf ; echo "##sample=$S" >> $OSS.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OSS.00.vcf
  bcftools norm -m - $OSS.vcf | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa | bedtools sort -header> $OSS.fix.vcf
 
  cat $OSS.fix.vcf | fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -rlen $HP_MTLEN -mfile $OS.max.vcf -ffile $OS.fix.vcf | \
    filterVcf.pl -sample $S -source $HP_M | bedtools sort  >> $OSS.00.vcf   #  to add -depth $HP_DP
  annotateVcf.sh $OSS.00.vcf
 
  #intersectVcf.pl $OS.00.vcf $OS.max.vcf | cat - $OSS.00.vcf |  uniqVcf.pl | bedtools sort -header > $OSS.00.vcf.tmp
  intersectVcf.pl $OS.00.vcf $OS.max.vcf | differenceVcf.pl - $OSS.00.vcf  | perl -ane 'if(/(.+):0\.\d+$/) { print "$1:1\n"} else { print }' | cat - $OSS.00.vcf | uniqVcf.pl | bedtools sort -header > $OSS.00.vcf.tmp
  mv $OSS.00.vcf.tmp $OSS.00.vcf
fi

rm -f $OS.bam*
rm -f $OSS.*idx $OSS.*tsv $OSS.*stat
