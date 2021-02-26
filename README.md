# HP : Heteroplasmy Pipeline # 


## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINIX 
    SOFTWARE PACKAGES: git, wget, tar, unzip, autoconf, gcc(v4.8.5+), java(v1.8.0+), perl(v5.16.3+) 
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa, samtools, bedtools, fastp, samblaster, bcftools, htslib, vcftools, st
    JAVA JARS:         gatk, mutserve, haplogrep
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/HP.git

### SETUP ENVIRONMENT ###
    
    $ export SDIR=$PWD/HP/scripts/   # set script directory variable
 
### INSTALL PIPELINE PREREQUISITES (optional) ###

    $ $SDIR/install_prerequisites.sh  
    
### CHECK INSTALL ###
  
    # if successfull => "Success message!"
    $ $SDIR/checkInstall.sh
    $ cat checkInstall.log

## USAGE ##

### COPY init.sh & GENERATE INPUT FILE  ###

    # go to the working directory

    # copy init.sh ; edit if necessary
    $ cp -i $SDIR/init.sh .
    $ nano init.sh

    # generate an imput file which contains the list of the BAM/CRAM files to be processed 
    $ ADIR=alignment directory path

    $ find $ADIR/  | ls2in.pl  | sort > in.txt

    $ head in.txt
      #sampleName    inputFile         outputPath/prefix
      sim.A          bams/sim.A.bam	out/sim.A/sim.A
      sim.B          bams/sim.B.bam	out/sim.B/sim.B
      sim.C          bams/sim.C.bam	out/sim.C/sim.C
      ...

### GENERATE INDEX AND COUNT FILES ###

     $ cut -f2 in.txt | sed "s|^|$SDIR/samtools.sh |" > samtools.all.sh
     ... or (MARCC)
     $ cut -f2 in.txt | sed "s|^|sbatch --partition=shared --time=4:0:0 $SDIR/samtools.sh |" > samtools.all.sh

     $ chmod a+x ./samtools.all.sh
     $ ./samtools.all.sh

### GENERATE PIPELINE SCRIPT ###

    $ $SDIR/run.sh in.txt out/  > filter.all.sh
    $ chmod u+x ./filter.all.sh

### RUN PIPELINE  ###

    $ ./filter.all.sh				
    ... or
    $ nohup ./filter.all.sh &	                               # run in the backgroud
    ... or (MARCC)
    $ sbatch --partition=shared --time=24:0:0 ./filter.all.sh    # run using a job scheduler

## OUTPUT ##

    TAB/VCF Files: 

### 1st ITTERATION ### 

    count.tab 
    {mutect2,mutserve}.{03,05,10}.{concat,merge[.sitesOnly]}.vcf       # SNVs; 3,5,10% minimum heteroplasmy thold
    {mutect2,mutserve}.{03,05,10}.tab                                  # SNV counts
    cont.tab                                                           # reads, mtDNA-CN counts

### 2nd ITTERATION ###

    mutect2.mutect2.{03,05,10}.{concat,merge[.sitesOnly]}.vcf	
    mutect2.mutect2.{03,05,10}.tab

## EXAMPLE 1 : Usage ##

    #script      inFile outDir   humanRef   MTRef   SNVcalller
    $ $SDIR/run.sh in.txt out/    
    ... or
    $ $SDIR/run.sh in.txt out/   hs38DH.fa  RSRS.fa mutserve


## EXAMPLE 2 : 3 simulated datasets  ##

### input file ###
    $ cat in.txt 
      bams/sim.A.bam
      bams/sim.B.bam
      bams/sim.C.bam 

### after running samtools.sh ###
    $ ls bams/*
      bams/sim.A.bam
      bams/sim.A.bam.bai
      bams/sim.A.idxstats
      bams/sim.A.count
     ...

### after running the pipeline ###

#### output directories : 1/sample ####
  
    $ ls out/
      sim.A/
      sim.B/
      sim.C/ 
      ...

#### vcf files : concat & merge ####

    $ cat mutect2.03.concat.vcf  
      ...
      #CHROM    POS  ID  REF  ALT  QUAL  FILTER                      INFO              FORMAT    SAMPLE    
      chrM      64   .   C    T    .     clustered_events;haplotype  SM=sim.A;HV       GT:DP:AF  0|1:67:1  
      chrM      73   .   A    G    .     PASS                        SM=sim.B;HV;NUMT  GT:DP:AF  0/1:52:1  
      chrM      73   .   A    G    .     PASS                        SM=sim.C;HV       GT:DP:AF  0/1:43:1  
      ...

    $ cat mutect2.03.merge.vcf  
      ...  
      #CHROM    POS  ID  REF  ALT  QUAL  FILTER  INFO               FORMAT    sim.A     chrM.B    chrM.C    
      chrM      64   .   C    T    .     .       AC=1;AN=2;HV       GT:DP:AF  0|1:67:1  .:.:.     .:.:.     
      chrM      73   .   A    G    .     .       AC=1;AN=2;HV;NUMT  GT:DP:AF  .:.:.     0/1:52:1  .:.:.     
      chrM      73   .   A    G    .     .       AC=2;AN=4;HV       GT:DP:AF  0|1:70:1  .:.:.     0/1:43:1  
      ...

    $ cat mutect2.03.merge.sitesOnly.vcf
      #CHROM    POS  ID  REF  ALT  QUAL  FILTER  INFO               
      chrM      64   .   C    T    .     .       AC=1;AN=2;HV       
      chrM      73   .   A    G    .     .       AC=1;AN=2;HV;NUMT  
      chrM      73   .   A    G    .     .       AC=2;AN=4;HV       
      ....
     
#### SNV counts ####

    $ cat mutect2.03.tab            
      Run                       H           h   S   s   I  i  Hp  hp  Sp  sp  Ip  ip  
      sim.A                     31          44  28  36  3  8  3   0   0   0   3   0   
      sim.B                     27          42  25  34  2  8  5   0   3   0   2   0   
      sim.C                     39          43  35  35  4  8  2   0   0   0   2   0   

#### Haplogroups ####

    $ cat mutect2.haplogroup.tab    
      Run                       haplogroup  
      sim.A                     A2+(64)     
      sim.B                     B2          
      sim.C                     C           

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
    S                             # homozygous SNPs, 3% heteroplasmy rate
    s                             # heterozygous SNPs
    I                             # homozygous INDELs
    i                             # heterozygous INDELs
    ...
    ?p                            # "p" suffix stands for non homopolimeric
