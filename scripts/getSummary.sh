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


#cat IN.txt  | perl -ane 'print "ls $ENV{D}/$F[1]/$F[1].\n";'

#read counts
cut -f3 $IN | sed "s|$|.count|" | xargs cat | uniq.pl -i 0  | sed 's|^sample|Run|' > $D/filter.tab

#cvg
cut -f3 in.txt | sed "s|$|.cvg.stat|" | xargs cat | uniq.pl -i 0  > $D/cvg.tab

#snv annotation
cut -f3 $IN | sed "s|$|.$M.00.vcf|" | xargs cat | bedtools sort -header | sed 's|rCRS|chrM|g'  > $D/$M.00.concat.vcf
cat $D/$M.00.concat.vcf | grep "^#" > $D/$M.00.concat.vcf.tmp
cat $D/$M.00.concat.vcf | grep -v "^#" | sort -k1,1 -k2,2n -k4,4 -k5,5 >>  $D/$M.00.concat.vcf.tmp
mv  $D/$M.00.concat.vcf.tmp  $D/$M.00.concat.vcf
annotateVcf.sh $D/$M.00.concat.vcf

#snv counts
snpCount.sh $IN $D $M $T1
snpCount.sh $IN $D $M $T2
snpCount.sh $IN $D $M $T3

#cleanup
rm -f fastp.html fastp.json

##########################################################

if [[ $M != "mutect2" ]] ; then exit 0 ; fi

#haplogroups
cut -f3 $IN | sed "s|$|.$M.haplogroup|" | xargs cat | grep -v SampleID | sed 's|"||g' | awk '{print $1,$2}' | sort -u | \
  perl -ane 'BEGIN {print "Run\thaplogroup\n"} print "$F[0]\t$F[1]\n";' | sed 's|\.MT||' | \
  tee $D/$M.haplogroup.tab | \
  perl -ane 'print "$1\n" if(/(^Run.+)/ or /(\S+\s+L\d)(.*)/ or /(\S+\s+HV)(.*)/ or /(\S+\s+JT)(.*)/ or /(\S+\s+\w)(.*)/);' > $D/$M.haplogroup1.tab

#fasta
cut -f3 $IN | sed "s|$|.$M.fa|"        | xargs cat > $D/$M.fa
samtools faidx  $D/$M.fa

#consensus cvg; tmp
#cut -f3 $IN | sed "s|$|.$M.merge.bed|"  | xargs cat > $D/$M.merge.bed
