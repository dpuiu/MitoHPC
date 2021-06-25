#!/usr/bin/bash -eux

##########################################################

#Program that summarizes aligned read counts
#Input arguments

IN=$1
D=$2                    # out dir
M=$3                    # snp caller: mutect2 or mutserve
T1=$4             	# Heteroplasmy tholds
T2=$5
T3=$6

##########################################################


cat IN.txt  | perl -ane 'print "ls $ENV{D}/$F[1]/$F[1].count\n";'

#read counts
if [ ! -s $D/filter.tab ] ; then
  find $D/*/ -name "*.count"   | sort | xargs cat | uniq.pl -i 0  | sed 's|^sample|Run|' > $D/filter.tab
fi 

#cvg
if [ ! -s $D/cvg.tab ] ; then
  find $D/*/ -name "*cvg.stat" | sort | xargs grep -v min -m 1 -H | egrep -v "mutect2|mutserve" | sort | perl -ane 'BEGIN { print join "\t",("Run","min","q1","q2","q3","max","mean\n") } print "$1\n" if(/.+:(.+)/);' > $D/cvg.tab
fi

#snv annotation
find $D/*/ -name "*.$M.00.vcf" -not -name "*.$M.$M.00.vcf" | sort | xargs cat | bedtools sort -header > $D/$M.00.concat.vcf
cat $D/$M.00.concat.vcf | grep "^#" > $D/$M.00.concat.vcf.tmp
cat $D/$M.00.concat.vcf | grep -v "^#" | sort -k1,1 -k2,2n -k4,4 -k5,5 >>  $D/$M.00.concat.vcf.tmp
mv  $D/$M.00.concat.vcf.tmp  $D/$M.00.concat.vcf
annotateVcf.sh $D/$M.00.concat.vcf

#snv counts
snpCount.sh $IN $D $M $T1
snpCount.sh $IN $D $M $T2
snpCount.sh $IN $D $M $T3

##########################################################

if [[ $M != "mutect2" ]] ; then exit 0 ; fi

#haplogroups
find $D/ -name "*.$M.haplogroup" | sort | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
  perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | sed 's|\.MT||' | \
  tee $D/$M.haplogroup.tab | \
  perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $D/$M.haplogroup1.tab

#fasta
find $D/*/ -name *.$M.fa        | sort | xargs cat > $D/$M.fa
samtools faidx  $D/$M.fa

#consensus cvg
find $D/*/ -name *.$M.merge.bed | sort | xargs cat > $D/$M.merge.bed

#cleanup
rm -f fastp.html fastp.json

