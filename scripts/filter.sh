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

I=$1
O=$2
M=$3
H=$4  ; test -s $H.fa ; test -s $H.NUMT.fa
FO=$5 ; test -s $FO.fa
F=$6  ; test -s $F.fa

########################################################################################################################################
#set variables

export N=`basename $O .mutect2` #2nd mutect2 itteration
export RO=`basename $FO .fa`
export R=`basename $F`
export R=`basename $R .mutect2`
export R=`basename $R .mutect2`

IDIR=`dirname $I`
ODIR=`dirname $O`; mkdir -p $ODIR
P=1                						# number of processors
MINSIZE=16500
MM=4g	# new (samtools : sorting mem mem)
MS=2g	# new (java vm gatk)
MX=2g
#########################################################################################################################################
#test input files
test -s $I
if [ ! -s $I.bai ] && [ ! -s $I.crai ] ; then exit 1 ; fi
if [ $(stat -c%s $F.fa) -lt $MINSIZE ] ; then exit 1 ; fi

#########################################################################################################################################
#format references

if [ ! -s $F.fa.fai ] ; then
  samtools faidx $F.fa
fi

if [ ! -s $F.dict   ] ; then
  java -jar $JDIR/gatk.jar CreateSequenceDictionary --REFERENCE $F.fa --OUTPUT $F.dict
fi

if [ ! -s $F+.fa  ] ; then
  cat $F.fa.fai | perl -ane 'print "$F[0]\t0\t$F[1]\n$F[0]\t0\t$ENV{E}\n";' | bedtools getfasta -fi $F.fa -bed - -fo /dev/stdout | grep -v "^>" | perl -ane 'BEGIN { print ">$ENV{R}\n" } ;print;' > $F+.fa
  cat $H.NUMT.fa >> $F+.fa
  cat $F+.fa | perl -ane 'if(/^>/) { print "\n" if($.>1); print} else { chomp ;print} END{print "\n"}' > $F+.norm.fa ; mv $F+.norm.fa $F+.fa
fi

if [ ! -s $F+.bwt ] ; then
  bwa index $F+.fa -p $F+
fi

#########################################################################################################################################
#filter & realign reads   ##   -F 20 (fwd); -f 0x10 (rev)
if  [ ! -s $O.bam.bai ] ; then
  test -s $I
  if [ ! -s $I.bai ] && [ ! -s $I.crai ] ; then exit 1 ; fi

  samtools view $I $MT $NUMT -bu -F 0x900 -T $H.fa | \
    samtools sort -n -O SAM -m $MM | \
    perl -ane 'if(/^@/) {print} elsif($P[0] eq $F[0]) {print $p,$_}; @P=@F; $p=$_;' | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout $FOPT -j $O.json  -h $O.html | \
    bwa mem $F+ - -p -v 1 -t $P -Y -R "@RG\tID:$N\tSM:$N\tPL:ILLUMINA" -v 1 | \
    samblaster --removeDups --addMateTags  | \
    perl -ane 'if(/^@/) { print } else { last if($P[0] ne $F[0] and $L>$ENV{L}) ;  print if($F[2] eq $ENV{R} and $F[6] eq "="); @P=@F ; $L++}' | \
    circSam.pl -ref_len $F.fa.fai | grep -v "^$" | \
    grep -v -P '^\@SQ\tSN:chr1' | \
    samtools view -bu | \
    samtools sort -m $MM > $O.bam

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
#get split alignments
samtools view -h $O.bam | sam2bedSA.pl | uniq.pl -i 3 | sort -k2,2n -k3,3n > $O.sa.bed

#########################################################################################################################################
#compute SNP/INDELs using mutect2 or mutserve

if [ ! -s $O.$M.vcf ] ; then
  if [ "$M" == "mutect2" ] ; then
    java -Xms$MS -Xmx$MX -jar $JDIR/gatk.jar Mutect2           -R $F.fa -I $O.bam                             -O $O.$M.vcf
    java -Xms$MS -Xmx$MX -jar $JDIR/gatk.jar FilterMutectCalls -R $F.fa -V $O.$M.vcf --min-reads-per-strand 2 -O $O.${M}F.vcf ; mv $O.${M}F.vcf  $O.$M.vcf ; rm $O.${M}F.vcf* $O.$M.vcf.*
    if [ -s $O.max.vcf ] ; then
      cat $O.max.vcf | fixsnpPos.pl -ref $RO -rfile $FO.fa -file /dev/stdin $O.$M.vcf > $O.${M}F.vcf ; mv $O.${M}F.vcf $O.$M.vcf
    fi
  elif [ "$M" == "mutserve" ] ; then
    if [ "$R" == "chrM" ] ||  [ "$R" == "rCRS" ] ||  [ "$R" == "RSRS" ] ; then
      java -Xms$MS -Xmx$MX -jar $JDIR/mutserve.jar call --deletions --insertions --level 0.01 --output $O.$M.vcf --reference $F.fa $O.bam
    else
      echo "Wrong mutserve reference"
      exit 1
    fi
  else
    echo "Unsuported SNV caller"
    exit 1
  fi
fi
test -s $O.$M.vcf

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

  #if [ $(stat -c%s " $O.$M.fa") -lt $MINSIZE ] ; then exit 1 ; fi

  cat $O.$M.fa | perl -ane 'if(/^>/) { print "\n" if($.>1); print} else { chomp ;print} END{print "\n"}' > $O.$M.norm.fa
  mv $O.$M.norm.fa $O.$M.fa
  bwa index $O.$M.fa  -p $O.$M
  bedtools bamtofastq -i $O.bam -fq /dev/stdout | bwa mem $O.$M - -v 1 -t $P -v 1 -k 63 | samtools sort -m $MM | bedtools bamtobed -tag NM -cigar | perl -ane 'print if($F[4]==0);' | bedtools merge -d -5  | bed2bed.pl > $O.$M.merge.bed
  rm -f $O.$M.*{sa,amb,ann,pac,bwt}

  ########################################################################################################################################
  # get haplogroup
  if [ "$R" == "rCRS" ] ||  [ "$R" == "chrM" ] ; then
    java -Xms$MS -Xmx$MX -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup
  elif [ "$R" == "RSRS" ] ; then
    java -Xms$MS -Xmx$MX -jar $JDIR/haplogrep.jar classify --in $O.$M.vcf  --format  vcf  --out $O.$M.haplogroup --rsrs
  fi
fi

if  [ -s $O.$M.fa ]  && [ ! -s $O.$M.haplocheck ] ; then
  java -Xms$MS -Xmx$MX -jar $JDIR/haplocheck.jar --out $O.$M.haplocheck $O.$M.vcf
fi
