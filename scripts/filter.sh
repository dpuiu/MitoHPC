#!/usr/bin/env bash
set -e

#######################################################################################################################################

#Program that runs the heteroplasmy pipeline on a single sample

#Input arguments
#1: sample names
#2: BAM/CRAM alignment file; full path
#3: output prefix; full path

########################################################################################################################################
#set variables

export S=$1
N=`basename $2 .bam`
export N=`basename $N .cram`
IDIR=`dirname $2`
I=$IDIR/$N
ODIR=`dirname $3`; mkdir -p $ODIR
O=$3

#########################################################################################################################################
#test input files

test -s $2

#test BAM/CRAM file sorted
samtools view -H $2 | grep -m 1 -P "^@HD.+coordinate$" > /dev/null

if [ ! -s $2.bai ] && [ ! -s $2.crai ]; then  samtools index -@ $HP_P $2 ; fi
if [ ! -s $I.idxstats ] ;               then  samtools idxstats -@ $HP_P $2 > $I.idxstats ; fi

if [ ! -s $I.count ] ; then
  cat $I.idxstats | idxstats2count.pl -sample $S -chrM $HP_RMT >  $I.count
  samtools view -F 0x900 $2 $HP_ENUMT -c -T $HP_RDIR/$HP_RNAME.fa | sed 's|^|NUMT\n|' | paste $I.count - > $I.count+ ; mv $I.count+ $I.count
fi

if [ $HP_I -lt 1 ] ; then exit 0 ; fi


#########################################################################################################################################
#sample reads

R=`tail -1 $I.count  | perl -ane '$S=$ENV{HP_L}/$F[3]; $S=0.9999 if($S>0.9999); print $S'`
if [ ! -s $O.fq ] ; then

  RCOUNT=`samtools view -H $2 | grep -c "^@SQ"`
  if [ $HP_RCOUNT != $RCOUNT ] ; then echo "ERROR: HP_RCOUNT=$HP_RCOUNT does not match the number of \@SQ lines=$RCOUNT in $2"; exit 1 ; fi

  samtools view -s $R $2 $HP_RMT $HP_RNUMT -bu -F 0x900 -T $HP_RDIR/$HP_RNAME.fa   | \
    samtools view -h | \
    samtools sort -n -O SAM -m $HP_MM | \
    perl -ane 'if(/^@/) {print} elsif($P[0] eq $F[0]) {print $p,$_}; @P=@F; $p=$_;' | \
    samblaster --removeDups --addMateTags | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout $HP_FOPT  > $O.fq
fi

#########################################################################################################################################
#realign reads

if  [ ! -s $O.bam ] ; then
  cat $O.fq | \
    bwa mem $HP_RDIR/${HP_MT}+ - -p -v 1 -t $HP_P -Y -R "@RG\tID:$1\tSM:$1\tPL:ILLUMINA" -v 1 | \
    circSam.pl -ref_len $HP_RDIR/$HP_MT.fa.fai | grep -v "^$" | \
    perl -ane 'next if(/^\@SQ\tSN:(\S+)/ and $1 ne $ENV{HP_MT}); next if(!/^\@/ and ($F[2] ne $ENV{HP_MT} or $F[6] ne "=")); print;' | \
    samtools view -bu | \
    samtools sort -m $HP_MM > $O.bam

  samtools index $O.bam
fi
#########################################################################################################################################
#count aligned reads; compute cvg; get stats; gets split alignments

if [ ! -s $O.count ]    ; then samtools idxstats $O.bam | idxstats2count.pl -sample $S -chrM $HP_MT > $O.count ; fi
if [ ! -s $O.cvg ]      ; then cat $O.bam | bedtools bamtobed -cigar | grep "^$HP_MT" | bedtools genomecov -i - -g $HP_RDIR/$HP_MT.fa.fai -d > $O.cvg ; fi
if [ ! -s $O.cvg.stat ] ; then cat $O.cvg | cut -f3 | st.pl  --summary --mean | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{S}\t$_" }'  > $O.cvg.stat ; fi
if [ ! -s $O.sa.bed ]   ; then samtools view -h $O.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed ; fi
#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

OM=$O.$HP_M
if [ ! -s $OM.00.vcf ] ; then
  if [ "$HP_M" == "mutect2" ] ; then
    java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $HP_RDIR/$HP_MT.fa -I $O.bam                            -O $OM.vcf-
    java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $HP_RDIR/$HP_MT.fa -V $OM.vcf- --min-reads-per-strand 2 -O $OM.vcf ; rm $OM.vcf?*
  elif [ "$HP_M" == "mutserve" ] ; then
    if [ "$HP_MT" == "chrM" ] ||  [ "$HP_MT" == "rCRS" ] ||  [ "$HP_MT" == "RSRS" ] ; then
      java $HP_JOPT -jar $HP_JDIR/mutserve.jar call --deletions --insertions --level 0.01 --output $OM.vcf --reference $HP_RDIR/$HP_MT.fa $O.bam
    else
      echo "Wrong mutserve reference"
      exit 1
    fi
  else
    echo "Unsuported SNV caller"
    exit 1
  fi

  test -s $OM.vcf

  # filter SNPs
  cat $HP_SDIR/$HP_M.vcf > $OM.00.vcf ; echo "##sample=$S" >> $OM.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OM.00.vcf
  cat $OM.vcf  | bcftools norm -m - | filterVcf.pl -sample $S -source $HP_M |  grep -v "^#" | sort -k2,2n -k4,4 -k5,5 | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa >> $OM.00.vcf
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

if  [ ! -s $O.fa.fai ]  ; then
  cat $OM.00.vcf | maxVcf.pl | bedtools sort -header |tee $OM.max.vcf | bgzip -f -c > $OM.max.vcf.gz ; tabix -f $OM.max.vcf.gz
  bcftools consensus -f $HP_RDIR/$HP_MT.fa $OM.max.vcf.gz | perl -ane 'chomp; if($.==1) { print ">$ENV{S}\n" } else { s/N//g; print } END {print "\n"}' > $OM.fa
  rm $OM.max.vcf.gz $OM.max.vcf.gz.tbi
  samtools faidx $OM.fa
fi

if [ $HP_I -lt 2 ] ; then exit 0 ; fi
if [ $HP_M != "mutect2" ] ; then exit 0 ; fi

########################################################################################################################################


#realign reads

if  [ ! -s $OM.bam.bai ] ; then
  samtools faidx $OM.fa $S:1-$HP_E | grep -v ">" | cat $OM.fa - > ${OM}+.fa
  cat $HP_RDIR/$HP_NUMT.fa >> ${OM}+.fa
  bwa index ${OM}+.fa -p ${OM}+

  cat $O.fq | \
    bwa mem ${OM}+ - -p -v 1 -t $HP_P -Y -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA" -v 1 | \
    circSam.pl -ref_len $OM.fa.fai | grep -v "^$" | \
    perl -ane 'next if(/^\@SQ\tSN:(\S+)/ and $1 ne $ENV{S}); next if(!/^\@/ and ($F[2] ne $ENV{S} or $F[6] ne "=")); print' | \
    samtools view -bu | \
    samtools sort -m $HP_MM > $OM.bam

  samtools index $OM.bam

  rm  ${OM}+.{sa,amb,ann,bwt,pac}
  #rm  ${OM}+*
fi

#########################################################################################################################################
#count aligned reads; compute cvg; get stats; gets split alignments

if [ ! -s $OM.count ]    ; then samtools idxstats $OM.bam | idxstats2count.pl -sample $S -chrM $S > $OM.count ; fi
if [ ! -s $OM.cvg ]      ; then cat $OM.bam | bedtools bamtobed -cigar | grep "^$HP_MT" | bedtools genomecov -i - -g $OM.fa.fai -d > $OM.cvg ; fi
if [ ! -s $OM.cvg.stat ] ; then cat $OM.cvg | cut -f3 | st.pl  --summary --mean | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{N}\t$_" }'  > $OM.cvg.stat ; fi
if [ ! -s $OM.sa.bed ]   ; then samtools view -h $OM.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $OM.sa.bed ; fi

#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

OMM=$OM.$HP_M

if [ ! -s $OMM.00.vcf ] ; then
  java $HP_JOPT -jar $HP_JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $OM.fa --OUTPUT $OM.dict

  java $HP_JOPT -jar $HP_JDIR/gatk.jar Mutect2           -R $OM.fa -I $OM.bam                            -O $OMM.vcf-
  java $HP_JOPT -jar $HP_JDIR/gatk.jar FilterMutectCalls -R $OM.fa -V $OMM.vcf- --min-reads-per-strand 2 -O $OMM.vcf+
  fixsnpPos.pl -ref $HP_MT -rfile $HP_RDIR/$HP_MT.fa -file $OM.max.vcf $OMM.vcf+ > $OMM.vcf ; rm $OMM.vcf?*

  test -s $OMM.vcf

  # filter SNPs
  cat $HP_SDIR/$HP_M.vcf > $OMM.00.vcf ; echo "##sample=$S" >> $OMM.00.vcf
  fa2Vcf.pl $HP_RDIR/$HP_MT.fa >> $OMM.00.vcf
  cat $OMM.vcf  | bcftools norm -m - | filterVcf.pl -sample $S -source $HP_M |  grep -v "^#" | sort -k2,2n -k4,4 -k5,5 | fix${HP_M}Vcf.pl -file $HP_RDIR/$HP_MT.fa >> $OMM.00.vcf
  annotateVcf.sh $OMM.00.vcf

  rm $OM.dict
fi

