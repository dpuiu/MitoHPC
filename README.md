# HP : Heteroplasmy Pipeline # 

# PREREQUISITES #

* bwa, samtools, bcftools, htslib, samblaster, vcftools
* picard, mutserve, gatk, haplogre
* hs38DH

# INSTALLATION #

## INSTALL PREREQUISITES ##

    wget https://netactuate.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
    wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
    wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2

    git clone git://github.com/GregoryFaust/samblaster.git
    git clone https://github.com/vcftools/vcftools.git

    wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
    wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
    wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip

## DOWNLOAD PIPELINE ##

    git clone https://github.com/dpuiu/HP.git

## SET PATHS & ENIRONMENTAL VARIABLES ##

    cd java/
    ln -s gatk-4.1.9.0.zip    gatk.jar
    ln -s haplogrep-2.2.9.jar haplogrepjar
    ln -s mutserve-1.3.4.jar  mutserve.jar  
    ln -s picard.jar        
    cd -

    cd scripts/
    export SDIR=`dirname $0`        # script directory
    export JDIR=$SDIR/../java/      # java jar directory
    export RDIR=$SDIR/../RefSeq/    # RefSeq directory ; should contain chrM.fa, rCRS.fa, RSRS.fa; link hs38DH.fa here
    export H=hs38DH.fa              # human reference
    export R=rCRS.fa                # or RSRS.fa
    export M=mutect2                # or mutserve
    export PATH=$SDIR:$PATH

    tree
    scripts
    |-- run.sh                    # main executable, calls "filter.sh" on multiple .bam/.cram files provided in an input file;
    |-- checkInstall.sh           # check prerequisites
    |-- init.sh                   # set environment variables
    |-- init_marcc.sh             # specific to MARCC
    |-- filter.sh                 # filter/realign reads, calls SNP/INDELs, filter SNP/INDELs at multiple heteroplamsy levels
    |-- readCount.sh              # count reads: all, mapped, chrM, filtered
    |-- snpCount.sh               # merge, count SNP/INDELs, HOM/HET(AF=) at multiple heteroplamsy levels
    |-- snpCount1.sh              # merge, count SNP/INDELs, HOM/HET(AF=) for a given heteroplamsy level
    |-- circSam.pl                # "circularizes" SAM alignments; extend reference, align & split reads spanning circ. point
    |-- count.pl                  # count values in a certain column (-i; 0 based)
    |-- fa2Vcf.pl                 # creates "##reference" & "##contig" VCF headers
    |-- filterVcf.pl              # filter VCF files; discard HETEROZYGOUS SNP/INDELs with AF less than a THOLD
    |-- fixmutect2Vcf.pl          # postprocess gatk Mutect2 output
    |-- fixmutserveVcf.pl         # postprocess mutserve output
    |-- join.pl                   # join 2 files by the 1st column
    |-- labelVcf.pl               # add the homopolimer tag(HP) to SNPs located at certain positions
    |-- maxVcf.pl                 # get the major allele
    |-- mutect2.vcf               # mutect2 VCF header
    |-- mutserve.vcf              # mutserve VCF header
    |-- uniq2.pl                  # filters unique lines based on 2 columns (-i 0 -j 1)
    java/                         # jars
    |-- gatk.jar
    |-- haplogrep.jar
    |-- mutserve.jar
    bin/                          # executables (in case they have not been already installed)
    |-- ...
    RefSeq/                       # references: chrM, hs38DH, rCRS
    |-- hs38DH.fa                
    |-- chrM.fa
    |-- rCRS.fa
    `-- RSRS.fa

## CHECK INSTALL ##

    checkInstall.sh

# USAGE #

### init (could be added to ~/.bash_profile) ###

    source HP/scripts/init.sh
    ... or (MARCC)
    source HP/scripts/init_marcc.sh

#### check install (once; if successfull => "Success message!") ####

    HP/scripts/checkInstall.sh

#### generate input file  ####

    find bams/ -name "*bam" > in.txt
    ... or
    find crams/ -name "*cram" > in.txt

#### split input file (optional; Ex: sets of 100) ####

    split -d -a 1 --numeric=1 -l 100 in.txt  in. --additional-suffix=.txt

#### generate pipeline script ####

    HP/scripts/run.sh in.txt   > filter.all.sh
    ... or
    HP/scripts/run.sh in.1.txt > filter.1.sh
    HP/scripts/run.sh in.2.txt > filter.2.sh
    HP/scripts/run.sh in.3.txt > filter.3.sh
    ...

#### run pipeline script  ####

    nohup ./filter.all.sh &
    ... or
    nohup ./filter.1.sh &
    nohup ./filter.2.sh &
    nohup ./filter.3.sh &
    ...
    ls filter.?.sh | parallel
    ... or (MARCC)
    sbatch --time=24:0:0 ./filter.all.sh

# OUTPUT #

### Files: ###
    mutect2.03.vcf
    mutect2.05.vcf
    mutect2.10.vcf
    mutect2.tab

# EXAMPLES #

#### use RSRS for realignment ####

    HP/scripts/run.sh in.txt . hs38DH.fa RSRS.fa

#### use rCRS for realignment, mutserve for SNP calling ####

    HP/scripts/run.sh in.txt . hs38DH.fa rCRS.fa mutserve

#### output ; simulated haplogroup A,B,C datasets ####

    head mutect2.03.vcf 
    ##fileformat=VCFv4.2
    ##source=Mutect2
    ##reference=file://../..//RefSeq//rCRS.fa>
    ##contig=<ID=rCRS,length=16569>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	                        INFO	                FORMAT	SAMPLE
    rCRS	64	.	C	T	.	clustered_events;haplotype	SNP;DP=69	        SM	chrM.A
    rCRS	73	.	A	G	.	clustered_events;haplotype	SNP;DP=72	        SM	chrM.A
    rCRS	73	.	A	G	.	PASS	                        SNP;DP=50	        SM	chrM.B
    rCRS	73	.	A	G	.	PASS	                        SNP;DP=43	        SM	chrM.C
    rCRS	146	.	T	C	.	clustered_events;haplotype	SNP;DP=88	        SM	chrM.A
    rCRS	153	.	A	G	.	clustered_events;haplotype	SNP;DP=89	        SM	chrM.A
    ...
    rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=71;HP	        SM	chrM.A
    rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=62;HP	        SM	chrM.B
    rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=57;HP	        SM	chrM.C
    ...
    rCRS	374	.	A	G	.	clustered_events	        SNP;DP=76;AF=0.139	SM	chrM.A
    rCRS	375	.	C	T	.	clustered_events;strand_bias	SNP;DP=65;AF=0.162	SM	chrM.B
    rCRS	378	.	C	T	.	clustered_events	        SNP;DP=81;AF=0.13	SM	chrM.C
    ...
 
    head mutect2.tab
    Run     haplogroup  03%S  03%s  03%I  03%i  05%S  05%s  05%I  05%i  10%S  10%s  10%I  10%i
    chrM.A  A2+(64)     28    35    3     8     28    35    3     8     28    34    3     8
    chrM.B  B2          25    34    2     8     25    34    2     8     25    34    2     8
    chrM.C  C           35    35    4     8     35    35    4     8     35    35    4     8

# LEGEND #

### metadata
    Run                           # SRR
    rdLen                         # AvgReadLength
    ...

### read counts

Run              
    all                           # all reads
    mapped                        # mapped reads	
    chrM                          # number of reads aligned to chrM
    filter                        # number of chrM used (subsample ~2000x cvg based on name)	

### computed coverage ####
    Gcvg                          # recomputed genome coverage: Bases/3217346917
    Mcvg                          # mitochondrion covearge: chrM*rdLen/16569

### mtDNA copy number ###
    M                             # 2*Mcvg/Gcvg

### mutect2 results ###
    haplogroup                    # mutect2 haploroup
    03%S                          # homozygous SNPs, 3% heteroplasmy rate
    03%S                          # heterozygous
    03%I                          # homozygous INDELs
    03%i                          # heterozygous INDELs
    ...
    05%
    10%

