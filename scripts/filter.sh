#!/bin/bash -eux

#######################################################################################################################################

#Program that runs the heteroplasmy pipeline on a single sample

#Input arguments
#1: I:  BAM/CRAM alignment file; full path
#2: O:  output prefix; full path
#3: M:  SNP calling method: mutect2(default) or mutserve
#4: H:  reference sequence path : hs38DH.fa assembly(default)
#5: FO: target sequence path : rCRS.fa assembly (default)    
#6: F:  target sequence path : rCRS.fa assembly (default) or sampleConsensus.fa(2nd itteration)

I=$1  ; test -s $I
O=$2
M=$3
H=$4  ; test -s $H
FO=$5 ; test -s $FO
F=$6  ; test -s $F

########################################################################################################################################
#set variables

export N=`basename $I .bam`        # sample name .bam or .cram
export N=`basename $N .cram`

export RO=`basename $FO .fa`
export R=`basename $F .fa`
export R=`basename $R .mutect2`
export R=`basename $R .mutect2`

G=${F%???} 

IDIR=`dirname $I`
ODIR=`dirname $O`; mkdir -p $ODIR
P=1                						# number of processors
MSIZE=16569

#########################################################################################################################################
#test input file

test -s $I
test -s $IDIR/$N.count
if [ ! -s $I.bai ] && [ ! -s $I.crai ] ; then exit 1 ; fi

if [ $(stat -c%s $F) -lt $MSIZE ] ; then exit 1 ; fi

#########################################################################################################################################
#format references

if [ ! -s $F.fai ]    ; then samtools faidx $F ; fi
if [ ! -s $G.dict   ] ; then java -jar $JDIR/picard.jar CreateSequenceDictionary R=$F O=$G.dict ; fi
if [ ! -s $G+.fa  ] ; then 
  cat $F.fai | perl -ane 'print "$F[0]\t0\t$F[1]\n$F[0]\t0\t$ENV{E}\n";' | bedtools getfasta -fi $F -bed - | grep -v "^>" | perl -ane 'BEGIN { print ">$ENV{R}\n" } ;print;' > $G+.fa
  cat $RDIR/$HG.NUMT.fa >> $G+.fa
  java -jar $JDIR/picard.jar NormalizeFasta I=$G+.fa O=$G+.norm.fa LINE_LENGTH=60 ;  mv $G+.norm.fa $G+.fa
fi
if [ ! -s $G+.bwt ] ; then bwa index $G+.fa -p $G+ ; fi

#########################################################################################################################################
#filter & realign reads

if  [ ! -s $O.bam ] ; then
  samtools view $I $MT $NUMT -bu -F 0x900 -T $H | \
    samtools sort -n -O SAM | \
    perl -ane 'if(/^@/) {print} elsif($L>$ENV{L}) {last} elsif($P[0] eq $F[0]) {print $p,$_ ; $L+=2}; @P=@F; $p=$_;' | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout | \
    bwa mem $G+ - -p -v 1 -t $P -Y -R "@RG\tID:$N\tSM:$N\tPL:ILLUMINA" -v 1 | \
    samblaster --removeDups --addMateTags  | \
    perl -ane 'print if(/^@/ or $F[2] eq $ENV{R} and $F[6] eq "=");' | \
    circSam.pl -ref_len $F.fai | grep -v "^$" | \
    grep -v -P '^\@SQ\tSN:chr1' | \
    samtools view -bu | \
    samtools sort > $O.bam
    samtools index $O.bam
fi
#########################################################################################################################################
#count aligned reads

if [ ! -s $O.count ] ; then
  cp $IDIR/$N.count $O.count
  samtools view $I $MT  -F 0x904 -c | awk '{print $1,"chrM" }'  >> $O.count
  samtools view $O.bam  -F 0x904 -c | awk '{print $1,"filter"}' >> $O.count
fi

#########################################################################################################################################
#get covearge at each chrM position ; get overall stats

if [ ! -s $O.cvg ] ; then
  cat $O.bam | bedtools bamtobed -cigar | bedtools genomecov -i - -g $F.fai -d > $O.cvg 
  cat $O.cvg | cut -f3 | st  --summary --mean | sed 's|^|'"$N"'\t|' > $O.cvg.stat
fi

#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

if [ ! -s $O.$M.vcf ] ; then
  if [ "$M" == "mutect2" ] ; then
    java -jar $JDIR/gatk.jar Mutect2           -R $F -I $O.bam                             -O $O.$M.vcf
    java -jar $JDIR/gatk.jar FilterMutectCalls -R $F -V $O.$M.vcf --min-reads-per-strand 2 -O $O.${M}F.vcf
    mv $O.${M}F.vcf  $O.$M.vcf ; rm $O.${M}F.vcf* $O.$M.vcf.*
    if [ -s $O.max.vcf ] ; then
      fixsnpPos.pl -ref $RO -rfile $FO -file $O.max.vcf $O.$M.vcf > $O.{$M}F.vcf
      mv $O.{$M}F.vcf $O.$M.vcf
    fi
  elif [ "$M" == "mutserve" ] && [ "$R" == "rCRS" ] ; then
    java -jar $JDIR/mutserve.jar analyse-local --input $O.bam --deletions  --insertions --level 0.01 --output $O.$M.vcf --reference $F
    cat $O.$M.vcf | perl -ane 'if(/^##/) { print } else { print join "\t",@F[0..9]; print "\n"}' | sed 's|^chrM|rCRS|g' | sed 's|.bam$||'  > $O.${M}F.vcf
    mv $O.${M}F.vcf $O.$M.vcf ; rm -f ${O}_raw.txt $O.txt
  fi

  ##########################################################################################################################################
  # filter SNPs
  cat $SDIR/$M.vcf > $O.$M.00.vcf
  fa2Vcf.pl $FO >> $O.$M.00.vcf
  cat $O.$M.vcf | bcftools norm -m -  | filterVcf.pl -sample $N -source $M | grep -v ^# | sort -k2,2n -k4,4 -k5,5 | fix${M}Vcf.pl -file $F   >> $O.$M.00.vcf
  #vcf-validator $O.$M.00.vcf
  cat $O.$M.00.vcf | filterVcf.pl -p 0.$T1  | tee $O.$M.$T1.vcf | filterVcf.pl -p 0.$T2  | tee $O.$M.$T2.vcf | filterVcf.pl -p 0.$T3  > $O.$M.$T3.vcf
fi

#########################################################################################################################################
#get new consensus

if  [ ! -s $O.fa ]  && [ ! -s $O.$M.fa ] ; then
  cat $O.$M.03.vcf | maxVcf.pl |  tee $O.$M.max.vcf  | bgzip -f -c > $O.$M.max.vcf.gz  ; tabix -f $O.$M.max.vcf.gz
  bcftools consensus -f $F $O.$M.max.vcf.gz | perl -ane 'if($.==1) { print ">$ENV{N}\n" } else { s/N//g; print }' > $O.$M.fa

  if [ $(stat -c%s " $O.$M.fa") -lt $MSIZE ] ; then exit 1 ; fi

  java -jar $JDIR/picard.jar NormalizeFasta I=$O.$M.fa O=$O.$M.norm.fa LINE_LENGTH=60 ; mv $O.$M.norm.fa $O.$M.fa ; rm $O.$M.max.vcf.*
  bwa index $O.$M.fa  -p $O.$M
  bedtools bamtofastq -i $O.bam -fq /dev/stdout | bwa mem $O.$M - -v 1 -t $P -v 1 -k 63 | samtools sort | bedtools bamtobed -tag NM -cigar | perl -ane 'print if($F[4]==0);' | bedtools merge -d -5  | bed2bed.pl > $O.$M.merge.bed
  rm -f $O.$M.*{sa,amb,ann,pac,bwt}

  ########################################################################################################################################
  # get haplogroup
  if [ "$R" == "rCRS" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup
  elif [ "$R" == "RSRS" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup --rsrs
  fi
fi

if [ -s $O.haplogroup ]; then
  cp $O.haplogroup $O.$M.haplogroup
fi
