#!/usr/bin/bash -eux

#######################################################################################################################################

#Program that runs the heteroplasmy pipeline on a single sample

#Input arguments
#1: I:  BAM/CRAM alignment file; full path
#2: O:  output prefix; full path
#3: M:  SNP calling method: mutect2(default) or mutserve
#4: H:  reference sequence path : hs38DH assembly(default)
#5: FO: target sequence path : rCRS (default)    
#6: F:  target sequence path : rCRS (default) or sampleConsensus(2nd itteration)

I=$1  ; test -s $I #tmp
O=$2
M=$3
H=$4  ; test -s $H.fa ; test -s $H.NUMT.fa
FO=$5 ; test -s $FO.fa
F=$6  ; test -s $F.fa


########################################################################################################################################
#set variables

export N=`basename $O .mutect2`        
export RO=`basename $FO .fa`
export R=`basename $F`
export R=`basename $R .mutect2`
export R=`basename $R .mutect2`

IDIR=`dirname $I`
ODIR=`dirname $O`; mkdir -p $ODIR
P=1                						# number of processors
MSIZE=16500

#########################################################################################################################################
#test input file

if [ ! -s $I.bai ] && [ ! -s $I.crai ] ; then exit 1 ; fi

if [ $(stat -c%s $F.fa) -lt $MSIZE ] ; then exit 1 ; fi

#########################################################################################################################################
#format references

if [ ! -s $F.fa.fai ] ; then
  samtools faidx $F.fa
fi

if [ ! -s $F.dict   ] ; then
  java -jar $JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $F.fa --OUTPUT $F.dict
fi

if [ ! -s $F+.fa  ] ; then
  cat $F.fa.fai | perl -ane 'print "$F[0]\t0\t$F[1]\n$F[0]\t0\t$ENV{E}\n";' | bedtools getfasta -fi $F.fa -bed - | grep -v "^>" | perl -ane 'BEGIN { print ">$ENV{R}\n" } ;print;' > $F+.fa
  cat $H.NUMT.fa >> $F+.fa
  #java -jar $JDIR/gatk.jar NormalizeFasta --INPUT $F+.fa --OUTPUT $F+.norm.fa --LINE_LENGTH 60
  cat $F+.fa | perl -ane 'if(/^>/) { print "\n" if($.>1); print} else { chomp ;print} END{print "\n"}' > $F+.norm.fa
  mv $F+.norm.fa $F+.fa
fi

if [ ! -s $F+.bwt ] ; then 
  bwa index $F+.fa -p $F+
fi

#########################################################################################################################################
#filter & realign reads

if  [ ! -s $O.bam.bai ] ; then
  samtools view $I $MT $NUMT -bu -F 0x900 -T $H.fa | \
    samtools sort -n -O SAM | \
    perl -ane 'if(/^@/) {print} elsif($L>$ENV{L}) {last} elsif($P[0] eq $F[0]) {print $p,$_ ; $L+=2}; @P=@F; $p=$_;' | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout | \
    bwa mem $F+ - -p -v 1 -t $P -Y -R "@RG\tID:$N\tSM:$N\tPL:ILLUMINA" -v 1 | \
    samblaster --removeDups --addMateTags  | \
    perl -ane 'print if(/^@/ or $F[2] eq $ENV{R} and $F[6] eq "=");' | \
    circSam.pl -ref_len $F.fa.fai | grep -v "^$" | \
    grep -v -P '^\@SQ\tSN:chr1' | \
    samtools view -bu | \
    samtools sort > $O.bam

  samtools index $O.bam
fi
#########################################################################################################################################
#count aligned reads

if [ ! -s $O.count ] ; then
  samtools idxstats $O.bam | idxstats2count.pl -sample $N -chrM $R > $O.count
fi

#########################################################################################################################################
#get covearge at each chrM position ; get overall stats

if [ ! -s $O.cvg.stat ] ; then
  cat $O.bam | bedtools bamtobed -cigar | bedtools genomecov -i - -g $F.fa.fai -d > $O.cvg
  cat $O.cvg | cut -f3 | st  --summary --mean | perl -ane 'if($.==1) { print "Run\t$_" } else { print "$ENV{N}\t$_" }'  > $O.cvg.stat
fi

#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

if [ ! -s $O.$M.vcf ] ; then
  if [ "$M" == "mutect2" ] ; then
    java -jar $JDIR/gatk.jar Mutect2           -R $F.fa -I $O.bam                             -O $O.$M.vcf
    java -jar $JDIR/gatk.jar FilterMutectCalls -R $F.fa -V $O.$M.vcf --min-reads-per-strand 2 -O $O.${M}F.vcf
    mv $O.${M}F.vcf  $O.$M.vcf ; rm $O.${M}F.vcf* $O.$M.vcf.*
    if [ -s $O.max.vcf ] ; then
      cat $O.max.vcf | fixsnpPos.pl -ref $RO -rfile $FO.fa -file /dev/stdin $O.$M.vcf > $O.${M}F.vcf
      mv $O.${M}F.vcf $O.$M.vcf
    fi
  elif [ "$M" == "mutserve" ] && [ "$R" == "rCRS" ] ; then
    java -jar $JDIR/mutserve.jar analyse-local --input $O.bam --deletions  --insertions --level 0.01 --output $O.$M.vcf --reference $F
    cat $O.$M.vcf | perl -ane 'if(/^##/) { print } else { print join "\t",@F[0..9]; print "\n"}' | sed 's|^chrM|rCRS|g' | sed 's|.bam$||'  > $O.${M}F.vcf
    mv $O.${M}F.vcf $O.$M.vcf ; rm -f ${O}_raw.txt $O.txt
  fi
fi
##########################################################################################################################################

if [ ! -s $O.$M.00.vcf ]; then
  # filter SNPs
  cat $SDIR/$M.vcf > $O.$M.00.vcf
  echo "##sample=$N" >> $O.$M.00.vcf
  fa2Vcf.pl $FO.fa >> $O.$M.00.vcf
  cat $O.$M.vcf  | bcftools norm -m - | filterVcf.pl -sample $N -source $M |  grep -v "^#" | sort -k2,2n -k4,4 -k5,5 | fix${M}Vcf.pl -file $F.fa >> $O.$M.00.vcf
fi
#########################################################################################################################################
#get new consensus

if  [ ! -s $O.fa ]  && [ ! -s $O.$M.fa ] ; then
  cat $O.$M.00.vcf | maxVcf.pl | bedtools sort -header |tee $O.$M.max.vcf | bgzip -f -c > $O.$M.max.vcf.gz ; tabix -f $O.$M.max.vcf.gz
  bcftools consensus -f $F.fa $O.$M.max.vcf.gz | perl -ane 'if($.==1) { print ">$ENV{N}\n" } else { s/N//g; print }' > $O.$M.fa
  rm $O.$M.max.vcf.gz $O.$M.max.vcf.gz.tbi

  #if [ $(stat -c%s " $O.$M.fa") -lt $MSIZE ] ; then exit 1 ; fi

  #java -jar $JDIR/gatk.jar NormalizeFasta --INPUT $O.$M.fa --OUTPUT $O.$M.norm.fa --LINE_LENGTH 60
  cat $O.$M.fa | perl -ane 'if(/^>/) { print "\n" if($.>1); print} else { chomp ;print} END{print "\n"}' > $O.$M.norm.fa
  mv $O.$M.norm.fa $O.$M.fa
  bwa index $O.$M.fa  -p $O.$M
  bedtools bamtofastq -i $O.bam -fq /dev/stdout | bwa mem $O.$M - -v 1 -t $P -v 1 -k 63 | samtools sort | bedtools bamtobed -tag NM -cigar | perl -ane 'print if($F[4]==0);' | bedtools merge -d -5  | bed2bed.pl > $O.$M.merge.bed
  rm -f $O.$M.*{sa,amb,ann,pac,bwt}

  ########################################################################################################################################
  # get haplogroup
  if [ "$R" == "rCRS" ] ||  [ "$R" == "chrM" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup
  elif [ "$R" == "RSRS" ] ; then
    java -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup --rsrs
  fi
fi


