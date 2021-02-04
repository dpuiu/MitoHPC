# HP : Heteroplasmy Pipeline # 

## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINIX 
    SOFTWARE PACKAGES: git, wget, tar, unzip, autoconf, gcc(v4.8.5+), java(v1.8.0+), perl(v5.16.3+) 
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa, samtools, bedtools, fastp, samblaster, bcftools, htslib, vcftools, st
    JAVA JARS:         picard, mutserve, gatk, haplogrep
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    git clone https://github.com/dpuiu/HP.git

### SETUP ENVIRONMENT ###
    
    export SDIR=$PWD/HP/scripts/   # set script directory variable
 
### INSTALL PIPELINE PREREQUISITES (optional) ###

    $SDIR/install_prerequisites.sh  
    
### CHECK INSTALL ###
  
    # if successfull => "Success message!"
    $SDIR/checkInstall.sh
    cat checkInstall.log

## USAGE ##

### COPY init.sh & GENERATE INPUT FILE  ###

    # go to the working directory

    # copy init.sh ; edit if necessary
    cp -i $SDIR/init.sh .
    nano init.sh

    # generate an imput file which contains the list of the BAM/CRAM files to be processed 
    # ADIR=alignment directory path
    find $ADIR -name "*.bam"  > in.txt
    ... or
    find $ADIR -name "*.cram" > in.txt
   
### GENERATE INDEX AND COUNT FILES ###

     sed "s|^|$SDIR/samtools.sh |" in.txt > samtools.all.sh
     ... or (MARCC)
     sed "s|^|sbatch --partition=shared --time=24:0:0 $SDIR/samtools.sh |" in.txt > samtools.all.sh

     chmod a+x ./samtools.all.sh
     ./samtools.all.sh

### GENERATE PIPELINE SCRIPT ###

    $SDIR/run.sh in.txt out/  > filter.all.sh
    chmod u+x ./filter.all.sh

### RUN PIPELINE  ###

    ./filter.all.sh				
    ... or
    nohup ./filter.all.sh &	                               # run in the backgroud
    ... or (MARCC)
    sbatch --partition=shared --time=24:0:0 ./filter.all.sh    # run using a job scheduler

## OUTPUT ##

    TAB/VCF Files: 

### 1st ITTERATION ### 

    count.tab 
 
    {mutect2,mutserve}.{03,05,10}.{concat,merge}.vcf.gz    # SNVs; 3,5,10% minimum heteroplasmy thold

    cont.tab                                               # reads, mtDNA-CN counts
    count_{mutect2,mutserve}.tab                           # reads, mtDNA-CN & SNV counts 

### 2nd ITTERATION ###

    mutect2.mutect2.{03,05,10}.concat.vcf.gz	
    count_mutect2.mutect2.tab

## EXAMPLE 1 : Usage ##

    #script      inFile outDir   humanRef   MTRef   SNVcalller
    $SDIR/run.sh in.txt out/    
    ... or
    $SDIR/run.sh in.txt out/     hs38DH.fa  RSRS.fa mutserve


## EXAMPLE 2 : 3 simulated datasets  ##

### input file ###
    cat in.txt 
      bams/sim.A.bam
      bams/sim.B.bam
      bams/sim.C.bam 

### after running samtools.sh ###
    ls bams/*
      bams/sim.A.bam
      bams/sim.A.bam.bai
      bams/sim.A.count
      bams/sim.B.bam
      bams/sim.B.bam.bai
      bams/sim.B.count
      bams/sim.C.bam
      bams/sim.C.bam.bai
      bams/sim.C.count
    
     cat bams/sim.A.count 
       11000 all
       11000 mapped
       10947 chrM

### after running the pipeline ###

    head count_mutect2.tab 
      Run    all        mapped     chrM    filter  M        haplogroup  03%S  03%s  03%I  03%i  05%S  05%s  05%I  05%i  10%S  10%s  10%I  10%i
      sim.A  740589366  739237125  487382  205095  256.05   A2+(64)     28    35    3     8     28    35    3     8     28    34    3     8
      sim.B  763658318  762297733  495743  205156  252.56   B2          25    34    2     8     25    34    2     8     25    34    2     8
      sim.C  749938586  748667901  590963  200121  306.55   C           35    35    4     8     35    35    4     8     35    35    4     8

    zcat mutect2.03.concat.vcf.gz 
      ##fileformat=VCFv4.2
      ##source=Mutect2
      ##reference=file://../..//RefSeq//rCRS.fa>
      ##contig=<ID=rCRS,length=16569>
      ##INFO=<ID=SM,Number=1,Type=String,Description="Sample">
      ##INFO=<ID=SNP,Number=0,Type=Flag,Description="SNP">
      ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="INDEL">
      ##INFO=<ID=HP,Number=0,Type=Flag,Description="Homoloplymer">
      ##INFO=<ID=HV,Number=0,Type=Flag,Description="Hypervariable">
      ...
      ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
      ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
      ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
      ...
      #CHROM  POS    ID  REF        ALT  QUAL  FILTER                        INFO                    FORMAT    SAMPLE
      rCRS    64     .   C           T    .     clustered_events;haplotype    SM=sim.A;SNP;HV         GT:DP:AF  0|1:67:1
      rCRS    73     .   A           G    .     clustered_events;haplotype    SM=sim.A;SNP;HV;NUMT    GT:DP:AF  0|1:70:1
      rCRS    73     .   A           G    .     PASS                          SM=sim.C;SNP;HV         GT:DP:AF  0/1:43:1
      rCRS    73     .   A           G    .     PASS                          SM=sim.B;SNP;HV         GT:DP:AF  0/1:52:1
      rCRS    146    .   T           C    .     clustered_events;haplotype    SM=sim.A;SNP;HV         GT:DP:AF  0|1:88:1
      rCRS    153    .   A           G    .     clustered_events;haplotype    SM=sim.A;SNP;HV         GT:DP:AF  0|1:89:1
      rCRS    235    .   A           G    .     clustered_events              SM=sim.A;SNP;HV         GT:DP:AF  0/1:74:1
      rCRS    247    .   GA          G    .     clustered_events              SM=sim.C;INDEL;HV       GT:DP:AF  0/1:67:1
      ...
      rCRS    285    .   CAA         C    .     clustered_events              SM=sim.C;INDEL;HV       GT:DP:AF  0/1:61:1
      rCRS    302    .   A           AC   .     clustered_events              SM=sim.A;INDEL;HV;HP    GT:DP:AF  0/1:71:1
      ...
      rCRS    375    .   C           T    .     clustered_events;strand_bias  SM=sim.B;SNP            GT:DP:AF  0/1:65:0.162
      rCRS    378    .   C           T    .     clustered_events              SM=sim.C;SNP            GT:DP:AF  0/1:81:0.13
      ...

    zcat mutect2.03.merge.vcf.gz
      ...
      ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">
      ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
      ##INFO=<ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">
      ...
      #CHROM  POS    ID  REF    ALT  QUAL  FILTER  INFO                     FORMAT    sim.A     sim.A           sim.C
      rCRS    64     .   C      T    .     .       AC=1;AN=2;SF=0f          GT:DP:AF  0|1:67:1   .              .
      rCRS    73     .   A      G    .     .       AC=3;AN=6;SF=0f,1,2      GT:DP:AF  0|1:70:1   0/1:52:1       0/1:43:1
      rCRS    146    .   T      C    .     .       AC=1;AN=2;SF=0f          GT:DP:AF  0|1:88:1   .              .
      rCRS    153    .   A      G    .     .       AC=1;AN=2;SF=0f          GT:DP:AF  0|1:89:1   .              .
      rCRS    235    .   A      G    .     .       AC=1;AN=2;SF=0f          GT:DP:AF  0/1:74:1   .              .
      rCRS    247    .   GA     G    .     .       AC=1;AN=2;INDEL;SF=2f    GT:DP:AF  .          .              0/1:67:1
      ...
      rCRS    285    .   CAA    C    .     .       AC=1;AN=2;INDEL;SF=2f    GT:DP:AF  .          .              0/1:61:1
      rCRS    302    .   A      AC   .     .       AC=1;AN=2;INDEL;SF=0f    GT:DP:AF  0/1:71:1   .              .
      ...
      rCRS    375    .  C       T    .     .       AC=1;AN=2;SF=1f          GT:DP:AF  .          0/1:65:0.162   .
      rCRS    378    .  C       T    .     .       AC=1;AN=2;SF=2f          GT:DP:AF  .          .              0/1:81:0.13
      ....

### LEGEND ###

    Run                           # Sample name

    all                           # number of reads in the sample reads 
    mapped                        # number of aligned reads in the sample
    chrM                          # number of reads aligned to chrM
    filter                        # number of chrM used (subsample ~2000x cvg based on name)

    Gcvg                          # recomputed genome coverage: Bases/3217346917
    Mcvg                          # mitochondrion covearge: chrM*rdLen/16569

    M                             # mtDNA copy number ; = 2*Mcvg/Gcvg

    haplogroup                    # sample haplogroup
    03%S                          # homozygous SNPs, 3% heteroplasmy rate
    03%s                          # heterozygous SNPs
    03%I                          # homozygous INDELs
    03%i                          # heterozygous INDELs
    ...
    05%
    10%
    ...
    ??%?p                         # "p" suffix stands for non homopimeric
