#!/bin/bash -eux

#1: I: input file.bam path/.bam
#2: O: output         path(ODIR)/prefix 
#3: M: input; snp callung method: mutect2 or mutserve
#4: H: input; human reference hs38DH.fa  path(RDIR)/.fa
#5: FO:input; original mito reference ; rCRS.fa    path(RDIR)/prefix.fa   
#6: F: input; original/new mito referece ; rCRS.fa  path(RDIR)/prefix.fa  path(ODIR)/prefix.fa 

I=$1  ; test -s $I
O=$2
M=$3
H=$4  ; test -f $H
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

ODIR=`dirname $O`; mkdir -p $ODIR

########################################################################################################################################
export HG=hs38DH
export MT=chrM
export NUMT='chr1:629084-634422 chr17:22521366-22521502 '   # chrM + 2 selected NUMT

#export HG=grch38_1kgmaj
#export MT=M                                                 # tmp; sim dataset
#export NUMT='1:629084-634422 17:22521366-22521502 '        


P=1                # number of processors
export L=222000    # ~2000x MT coverage
export E=300       # extension(circularization)

#########################################################################################################################################
#test input file

test -s $I
test -s $I.count
if [ ! -s $I.bai ] && [ ! -s $I.crai ] ; then exit 1 ; fi

#########################################################################################################################################
#format references

if [ ! -s $F.fai ]    ; then samtools faidx $F ; fi
if [ ! -s $G.dict   ] ; then java -jar $JDIR/picard.jar CreateSequenceDictionary R=$F O=$G.dict ; fi
if [ ! -s $G+.fa  ] ; then 
  cat $F.fai | perl -ane 'print "$F[0]\t0\t$F[1]\n$F[0]\t0\t$ENV{E}\n";' | bedtools getfasta -fi $F -bed - | grep -v "^>" | perl -ane 'BEGIN { print ">$ENV{R}\n" } ;print;' > $G+.fa
  cat $RDIR/$HG.NUMT.fa >> $G+.fa
  java -jar $JDIR/picard.jar NormalizeFasta I=$G+.fa O=$G+.fa2 LINE_LENGTH=60 ;  mv $G+.fa2 $G+.fa
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
  cp $I.count $O.count
  samtools view $I $MT  -F 0x900 -c | awk '{print $1,"chrM" }'  >> $O.count
  samtools view $O.bam  -F 0x900 -c | awk '{print $1,"filter"}' >> $O.count
fi

#########################################################################################################################################
#get covearge at each chrM position ; get overall stats

if [ ! -s $O.cvg ] ; then
  cat $O.bam | bedtools bamtobed -cigar | bedtools genomecov -i - -g $F.fai -d > $O.cvg 
  cat $O.cvg | cut -f3 | st  --summary  | sed 's|^|'"$N"'\t|' > $O.cvg.stat
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
  fa2Vcf.pl $F >> $O.$M.00.vcf
  cat $O.$M.vcf | bcftools norm -m -  | filterVcf.pl -sample $N -source $M | grep -v ^# | sort -k2,2n -k4,4 -k5,5 | fix${M}Vcf.pl -file $F   >> $O.$M.00.vcf
  cat $O.$M.00.vcf | bcftools norm -m +  | perl -ane 'if(/^#/) {print} else { chomp; @F=split /\t/; $F[6]=~s/multiallelic// unless($F[4]=~/,/); $F[6]="." unless($F[6]); print join "\t",@F;print "\n"}' | bcftools norm -m - > $O.$M.00.vcf2
  mv $O.$M.00.vcf2 $O.$M.00.vcf ; vcf-validator $O.$M.00.vcf
  cat $O.$M.00.vcf | filterVcf.pl -p 0.03  | tee $O.$M.03.vcf | filterVcf.pl -p 0.05  | tee $O.$M.05.vcf | filterVcf.pl -p 0.10  > $O.$M.10.vcf
fi

#########################################################################################################################################
#get new consensus

if  [ ! -s $O.fa ]  && [ ! -s $O.$M.fa ] ; then
  cat $O.$M.03.vcf | maxVcf.pl |  tee $O.$M.max.vcf  | bgzip -f -c > $O.$M.max.vcf.gz  ; tabix -f $O.$M.max.vcf.gz
  bcftools consensus -f $F $O.$M.max.vcf.gz | perl -ane 'if($.==1) { print ">$ENV{N}\n" } else { s/N//g; print }' > $O.$M.fa
  java -jar $JDIR/picard.jar NormalizeFasta I=$O.$M.fa O=$O.$M.fa2 LINE_LENGTH=60 ; mv $O.$M.fa2 $O.$M.fa ; rm $O.$M.max.vcf.*
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
