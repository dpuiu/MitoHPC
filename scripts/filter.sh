#!/bin/bash -eux

#1: input file.bam path/prefix.bam
#2: input hg19.fa  path/prefix.fa
#3: input rCRS.fa  path/prefix.fa
#4: output         path/prefix
#5: method: mutect2 or mutserve

I=$1 ; test -s $I
H=$2 ; test -s $H
F=$3 ; test -s $F
D=$4 ; mkdir -p $D
M=$5

########################################################################################################################################
#set variables

export N=`basename $I .bam` ; export N=`basename $N .cram`   # sample name .bam or .cram
export R=`basename $F .fa`                                   # rCSRC or RSRS or ...
export R=`basename $R .mutect2`
O=$D/$N                                                      # output dir.sample name
G=${F%???}                                                   # ...

export HG=hs38DH
export MT=chrM
export NUMT='chr1:629084-634422 chr17:22521366-22521502 '   # chrM + 2 selected NUMT

P=1                # number of processors
export L=222000    # ~2000x MT coverage
#export L=11100    # ~100x 
export E=255       # extension(circularization) ; 149 for 150bp reads

#########################################################################################################################################
#format references

if [ ! -s $I.bai ]    ; then samtools index $I ; fi
if [ ! -s $F.fai ]    ; then samtools faidx $F ; fi
if [ ! -s $G.dict   ] ; then java -jar $JDIR/picard.jar CreateSequenceDictionary R=$F O=$G.dict ; fi
if [ ! -s $G+$E.fa  ] ; then 
  cat $F.fai | perl -ane 'print "$F[0]\t0\t$F[1]\n$F[0]\t0\t$ENV{E}\n";' | bedtools  getfasta -fi $F -bed - | grep -v "^>" | perl -ane 'BEGIN { print ">$ENV{R}\n" } ;print;' > $G+$E.fa
  java -jar $JDIR/picard.jar NormalizeFasta I=$G+$E.fa O=$G+$E.fa2 LINE_LENGTH=60
  mv $G+$E.fa2 $G+$E.fa
fi
if [ ! -s $G+$E.bwt ] ; then bwa index $G+$E.fa -p $G+$E ; fi

#########################################################################################################################################
#filter & realign reads

if  [ ! -s $O.bam ] ; then
  samtools view $I $MT $NUMT -bu -F 0x900 -T $H | \
    samtools sort -n -O SAM | \
    perl -ane 'print if(/^@/ or $F[2] eq $ENV{MT} or $F[6] eq $ENV{MT});' | \
    perl -ane 'if(/^@/) {print} elsif($L>$ENV{L}) {last} elsif($P[0] eq $F[0]) {print $p,$_ ; $L+=2}; @P=@F; $p=$_;' | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout | \
    bwa mem $G+$E - -p -v 1 -t $P -Y -R "@RG\tID:$N\tSM:$N\tPL:ILLUMINA" -v 1 | \
    samblaster --removeDups --addMateTags  | \
    circSam.pl -ref_len $F.fai | grep -v "^$" | tee $O.sam | \
    samtools view -bu | \
    samtools sort > $O.bam
    samtools index $O.bam
fi

#########################################################################################################################################
#count aligned raeds

if [ ! -s $O.count ] ; then
  samtools view $I -c               | awk '{print $1,"all"}'    >  $O.count
  samtools view $I -F 0x904 -c      | awk '{print $1,"mapped"}' >> $O.count
  samtools view $I $MT  -F 0x900 -c | awk '{print $1,"chrM" }'  >>  $O.count
  samtools view $O.bam  -F 0x900 -c | awk '{print $1,"filter"}' >> $O.count
fi

#########################################################################################################################################
#get covearge at each chrM position ; get overall stats

if [ ! -s $O.cvg ] ; then
  cat $O.bam | \
    bedtools bamtobed -cigar | tee $O.bed | \
    bedtools genomecov -i - -g $F.fai -d > $O.cvg 
  cat $O.cvg | cut -f3 | st  --summary  | sed 's|^|'"$N"'\t|' > $O.cvg.stat
fi

#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

if [ ! -s $O.$M.vcf ] ; then
  if [  "$M" == "mutect2" ] ; then
    java -jar $JDIR/gatk.jar Mutect2           -R $F -I $O.bam                             -O $O.$M.vcf
    java -jar $JDIR/gatk.jar FilterMutectCalls -R $F -V $O.$M.vcf --min-reads-per-strand 2 -O $O.${M}F.vcf # --min-allele-fraction 0.03 --unique 2
    mv $O.${M}F.vcf  $O.$M.vcf
    rm $O.${M}F.vcf* $O.$M.vcf.*
  elif [ "$M" == "mutserve" ] ; then
    if [ "$R" == "rCRS" ] ; then
      java -jar $JDIR/mutserve.jar analyse-local --input $O.bam --deletions  --insertions --level 0.01 --output $O.$M.vcf --reference $F
      cat $O.$M.vcf | perl -ane 'if(/^##/) { print } else { print join "\t",@F[0..9]; print "\n"}' | sed 's|.bam$||'  > $O.${M}F.vcf
      mv $O.${M}F.vcf $O.$M.vcf
      rm -f ${O}_raw.txt $O.txt
    fi
  fi
fi
if [ ! -s $O.$M.vcf ] ; then exit 1 ; fi

#########################################################################################################################################
#filter SNP/INDELs at 1,3,5,10% heteroplasmy
##slippage|weak_evidence|strand_bias|germline|strict_strand|base_qual  ... germline|fragment|position

if  [ ! -s $O.$M.10.vcf ] ; then
  cat $SDIR/$M.vcf > $O.$M.00.vcf
  fa2Vcf.pl $F >> $O.$M.00.vcf
  cat $O.$M.vcf | bcftools norm -m -  | egrep -v 'strict_strand|weak_evidence|base_qual' | filterVcf.pl -sample $N -source $M | grep -v ^# | sort -k2,2n -k4,4 -k5,5 | fix${M}Vcf.pl -file $F | labelVcf.pl  >> $O.$M.00.vcf
  cat $O.$M.00.vcf | bcftools norm -m +  | perl -ane 'if(/^#/) {print} else { chomp; @F=split /\t/; $F[6]=~s/multiallelic// unless($F[4]=~/,/); $F[6]="." unless($F[6]); print join "\t",@F;print "\n"}' | bcftools norm -m - > $O.$M.00.vcf2
  mv $O.$M.00.vcf2 $O.$M.00.vcf
  vcf-validator $O.$M.00.vcf
  cat $O.$M.00.vcf | filterVcf.pl -p 0.03  | tee $O.$M.03.vcf | filterVcf.pl -p 0.05  | tee $O.$M.05.vcf | filterVcf.pl -p 0.10  > $O.$M.10.vcf
fi

#########################################################################################################################################
#get new consensus

if  [ ! -s $O.$M.fa ] ; then
  cat $O.$M.10.vcf | maxVcf.pl | tee $O.$M.max.vcf | bgzip -f -c > $O.$M.max.vcf.gz  ; tabix -f $O.$M.max.vcf.gz
  bcftools consensus -f $F $O.$M.max.vcf.gz | perl -ane 'if($.==1) { print ">$ENV{N}\n" } else { s/N//g; print }' > $O.$M.fa
  java -jar $JDIR/picard.jar NormalizeFasta I=$O.$M.fa O=$O.$M.fa2 LINE_LENGTH=60 ; mv $O.$M.fa2 $O.$M.fa
  rm $O.$M.max.vcf.*
fi

#########################################################################################################################################
#get haplogroup

if  [ ! -s $O.$M.haplogroup ] ; then
  if [ "$R" == "rCRS" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup
  elif [ "$R" == "RSRS" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup --rsrs
  fi
fi


