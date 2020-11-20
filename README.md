# PREREQUISITES #

### executables: bwa, samtools, bcftools, htslib, samblaster, vcftools ###
    wget https://netactuate.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
    wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
    wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
    git clone git://github.com/GregoryFaust/samblaster.git
    git clone https://github.com/vcftools/vcftools.git

### java: picard, mutserve, gatk, haplogrep ###
    wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
    wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
    wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip

### RefSeq: hs38DH ###
    wget ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens/hs38DH/hs38DH.fa

# FILES #

    $ tree <br>
    scripts/<br>
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
    bin/                          # executables (in case they have not been already installed)<br>
    |-- ...
    RefSeq/                       # references: chrM, hs38DH, rCRS
    |-- hs38DH.fa                 # to be downloaded separately
    |-- chrM.fa
    |-- rCRS.fa
    `-- RSRS.fa

# LEGEND #

### metadata ###
    Run                           # SRR
    rdLen                         # AvgReadLength
    ...

### read counts ###
    chrM                          # number of reads aligned to chrM

### computed coverage ####
    Gcvg                          # recomputed genome coverage: Bases/3217346917
    Mcvg                          # mitochondrion covearge: chrM*rdLen/16569

### mtDNA copy number ###
    M                             # Gcvg based:  2*Mcvg/Gcvg

### mutect2 results ###
    haplogroup                    # mutect2 haploroup
    03%S                          # homozygous SNPs, 3% heteroplasmy rate
    03%S                          # heterozygous
    03%I                          # homozygous INDELs
    03%i                          # heterozygous INDELs
    ...
    05%
    10%

# EXAMPLE 1 #

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
    HP/scripts/run.sh in.txt out > filter.all.sh
    ... or
    HP/scripts/run.sh in.1.txt out/ > filter.1.sh
    HP/scripts/run.sh in.2.txt out/ > filter.2.sh
    HP/scripts/run.sh in.3.txt out/ > filter.3.sh
    ...

#### run pipeline  ####
    nohup ./filter.all.sh &
    ... or
    nohup ./filter.1.sh &
    nohup ./filter.2.sh &
    nohup ./filter.3.sh &
    ...
    ls filter.?.sh | parallel
    ... or (MARCC)
    sbatch --time=24:0:0 ./filter.all.sh

# EXAMPLE 2 #

#### use RSRS for realignment ####
    HP/scripts/run.sh filter.all.txt filter.all/ hs38DH.fa RSRS.fa

# EXAMPLE 3 #
#### use rCRS for realignment, mutserve for SNP calling ####
    HP/scripts/run.sh filter.all.txt filter.all/ hs38DH.fa rCRS.fa mutserve
