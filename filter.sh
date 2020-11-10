#!/bin/bash -eux
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH --time=1:0:0
#SBATCH --ntasks-per-node=1

#1: input file; .bam
#2: MT referece

test -s $1
test -s $1.bai
N=`basename $1 .bam`   # sample name
Q=${N}F	       # output prefix
#Q=$N

M=$2
R=hs38DH
F='chr1:629084-634422 chr17:22521366-22521502 chrM'   # chrM + 2 selected NUMT

P=1               # number of processors
export L=222000   # ~2000x MT coverage
export L=11100    # ~100x 
E=255             # extension(circularization) ; 149 for 150bp reads

########################################################################################################################################
#test input files exist

HDIR=$(dirname $(readlink -f $0))/../
RDIR=$HDIR/RefSeq/
JDIR=$HDIR/java/

test -f $RDIR/$R.fa                                 # needed for unfiltered .cram input
test -f $RDIR/bwa/${M}+${E}.bwt
test -s $RDIR/$M.fa
test -s $RDIR/$M.dict
########################################################################################################################################

#test executables and Java jars
which perl		#usually available on Linux
which gcc		#MARCC:"module load gcc"
which java		#MARCC:"module load java
which bwa		#MARCC:"module load bwa"  or install from "https://sourceforge.net/projects/bio-bwa/files/"
which samtools		#install from "http://www.htslib.org/download/"
which bedtools		#MARCC:"module load bedtools"
which fastp		#install from "https://github.com/OpenGene/fastp"
which samblaster	#install from "https://github.com/GregoryFaust/samblaster"
which circSam.pl	#available under scripts
which filterVcf.pl	#available under scripts
which st		#simple statistics: install from "https://github.com/nferraz/st"

test -f $JDIR/gatk.jar
test -f $JDIR/haplogrep.jar
test -f $JDIR/mutserve.jar
#########################################################################################################################################

if  [ ! -s $N.bam ] ; then
  samtools view $1 $F -b -F 0x100 -F 0x800 -T $RDIR/$R.fa > $N.bam
  samtools index $N.bam
  samtools view $N.bam | cut -f3 | uniq -c > $N.count
fi

#########################################################################################################################################
#filter alignments

if [ ! -s $Q.bam ] ; then
  samtools view $N.bam -bu | \
    samtools sort -n -O SAM | \
    perl -ane 'print if(/^@/ or $F[2] eq "chrM" or $F[6] eq "chrM");' | \
    perl -ane 'if(/^@/) {print} elsif($L>$ENV{L}) {last} elsif($P[0] eq $F[0]) {print $p,$_ ; $L+=2}; @P=@F; $p=$_;' | \
    samtools view -bu | \
    bedtools bamtofastq  -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout | \
    fastp --stdin --interleaved_in --stdout | \
    bwa mem ${RDIR}/bwa/${M}+${E} - -p -v 1 -t $P -Y -R "@RG\tID:$Q\tSM:$Q\tPL:ILLUMINA" -v 1 | \
    samblaster --removeDups --addMateTags  | \
    circSam.pl -ref_len $RDIR/$M.fa.fai | \
    samtools view -bu | \
    samtools sort > $Q.bam
    samtools index $Q.bam
fi

#########################################################################################################################################
#need to install "st"
#get cvg at each chrM position ; overall stats
if [ ! -s $Q.cvg ] ; then
  cat $Q.bam | \
    bedtools bamtobed -cigar | tee $Q.bed | \
    bedtools genomecov -i - -g $RDIR/$M.fa.fai -d > $Q.cvg 
    cat $Q.cvg | cut -f3 | st  --summary  | sed 's|^|'"$Q"'\t|' > $Q.cvg.stat
fi

#########################################################################################################################################
#compute SNPs, MT consensus using mutect2/mutserve ; filter 3% heteroplasmy 

if [ ! -s $Q.mutect2.vcf ] ; then
  java -jar $JDIR/gatk.jar Mutect2 -R $RDIR/$M.fa -I $Q.bam  -O $Q.mutect2.vcf
  java -jar $JDIR/gatk.jar FilterMutectCalls -R $RDIR/$M.fa  -V $Q.mutect2.vcf --min-reads-per-strand 2  -O $Q.mutect2F.vcf 
  mv $Q.mutect2F.vcf $Q.mutect2.vcf
  rm $Q.mutect2F.vcf* $Q.mutect2.vcf.idx $Q.mutect2.vcf.stats
  annotateVcf.sh $Q mutect2 $M
fi
#########################################################################################################################################

if [ "$M" == "rCRS" ] || [ "$M" == "RSRS" ] ; then
  if [ ! -s $Q.mutserve.vcf ] ; then
    java -jar $JDIR/mutserve.jar analyse-local --input $Q.bam --deletions  --insertions --level 0.01 --output $Q.mutserve.vcf --reference $RDIR/$M.fa
    rm -r ${Q}_raw.txt $Q.txt
    annotateVcf.sh $Q mutserve $M
  fi
fi
