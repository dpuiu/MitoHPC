# HP : Heteroplasmy Pipeline # 


## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINIX 
    SOFTWARE PACKAGES: git, wget, tar, unzip, autoconf, gcc(v4.8.5+), java(v1.8.0+), perl, python 
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa, samtools, bedtools, fastp, samblaster, bcftools, htslib, vcftools, st
    JAVA JARS:         gatk, mutserve, haplogrep, haplocheck
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/HP.git

### SETUP ENVIRONMENT ###
    
    $ PWD=`pwd`
    $ export HP_SDIR=path_to_pipeline_home_directory/HP/scripts/   # set script directory variable

### INSTALL PIPELINE PREREQUISITES (if necessary) ###

    $ sudo $HP_SDIR/install_sysprerequisites.sh
    $ $HP_SDIR/install_prerequisites.sh  
    
### CHECK INSTALL ###
  
    # if successfull => "Success message!"

    $ $HP_SDIR/checkInstall.sh
    $ cat checkInstall.log

## USAGE ##

### COPY init.sh SCRIPT ###

    # go to the working directory ; copy init.sh 

    $ cp -i $HP_SDIR/init.sh .

### SETUP ENVIRONMENT ###

    # set variables (in terminal or init.sh)

    $ nano init.sh 
        ...
       export HP_ADIR=$PWD/bams/                                           # alignment dir
       export HP_ODIR=$PWD/out/ ;                                          # output dir  
       export HP_IN=$PWD/in.txt                                            # input file
    
    $ . ./init.sh

    $ mkdir -p $HP_ODIR                                       

    $ find $HP_ADIR/ -type f -name "*.bam" -o -name "*.cram" | \
        $HP_SDIR/ls2in.pl -out $HP_ODIR | sort > $HP_IN

### GENERATE ALIGNMENT INDEX AND READ COUNT FILES ###

     $ cut -f2 $HP_IN  | sed "s|^|$SH $HP_SDIR/samtools.sh |" > samtools.all.sh
     $ . ./samtools.all.sh

### COMPUTE mtDNA-CN ####
    
     $ cut -f2 $HP_IN | sed -r 's|(.*)\.|\1\t|g' | cut -f1 | sed 's|$|.count|' | xargs cat | \
         $HP_SDIR/uniq.pl | $HP_SDIR/getCN.pl > $HP_ODIR/count.tab
  
### RUN PIPELINE  ###
 
    $ $HP_SDIR/run.sh > filter.all.sh
    $ . ./filter.all.sh                                             
     
## OUTPUT ##

    under $ODIR:

    VCF/TAB/SUMMARY/FASTA Files: 

### 1st ITTERATION ### 

    {mutect2,mutserve}.{03,05,10}.{concat,merge[.sitesOnly]}.vcf     # SNVs; 3,5,10% heteroplasmy thold
    {mutect2,mutserve}.{03,05,10}.tab                                # SNV counts
    {mutect2,mutserve}.{03,05,10}.summary                            # SNV count summaries

    {mutect2,mutserve}.fa                                            # consensus sequence
    {mutect2,mutserve}.haplogroup[1].tab                             # haplogroup
    {mutect2,mutserve}.haplocheck.tab                                # contamination screen   
    count.tab                                                        # reads  & mtDNA-CN counts
    cvg.tab                                                          # coverage stats
    

### 2nd ITTERATION ###

    mutect2.mutect2.{03,05,10}.{concat,merge[.sitesOnly]}.vcf	
    mutect2.mutect2.{03,05,10}.tab

## noNUMT ##

   removed the NUMT labeled SNV's

## EXAMPLE 1 : Usage ##

    $ $HP_SDIR/run.sh     

## EXAMPLE 2 : 3 simulated datasets  ##

### input file ###

    $ head $IN
      #sampleName     inputFile          outputPath/prefix
      chrM.A          bams/chrM.A.bam    out/chrM.A/chrM.A
      chrM.B          bams/chrM.B.bam    out/chrM.B/chrM.B
      chrM.C          bams/chrM.C.bam    out/chrM.C/chrM.C
      ...

 ### read counts and mtDNA-CN(M)

     $ head $ODIR/count.tab
       Run	 all        mapped     chrM    M
       chrM.A    851537886  848029490  396766  181.7
       chrM.B    884383716  882213718  506597  223.01
       chrM.C    786560467  785208588  503241  248.9
       ...


### after running samtools.sh ###
    $ ls $ADIR/*
      bams/chrM.A.bam
      bams/chrM.A.bam.bai
      bams/chrM.A.idxstats
      bams/chrM.A.count
     ...

### after running the pipeline ###

#### output directories : 1/sample ####
  
    $ ls out/
      chrM.A/
      chrM.B/
      chrM.C/ 
      ...

#### vcf files : concat & merge ####

    $ cat mutect2.03.concat.vcf  
      ...
      #CHROM  POS ID  REF  ALT  QUAL FILTER                      INFO               FORMAT    SAMPLE    
      chrM    64  .   C    T    .    clustered_events;haplotype  SM=chrM.A;HV       GT:DP:AF  0|1:67:1  
      chrM    73  .   A    G    .    PASS                        SM=chrM.B;HV;NUMT  GT:DP:AF  0/1:52:1  
      chrM    73  .   A    G    .    PASS                        SM=chrM.C;HV       GT:DP:AF  0/1:43:1  
      ...

    $ cat mutect2.03.merge.vcf  
      ...  
      #CHROM  POS ID  REF  ALT  QUAL FILTER  INFO               FORMAT    chrM.A     chrM.B    chrM.C    
      chrM    64  .   C    T    .    .       AC=1;AN=2;HV       GT:DP:AF  0|1:67:1   .:.:.     .:.:.     
      chrM    73  .   A    G    .    .       AC=1;AN=2;HV;NUMT  GT:DP:AF  .:.:.      0/1:52:1  .:.:.     
      chrM    73  .   A    G    .    .       AC=2;AN=4;HV       GT:DP:AF  0|1:70:1   .:.:.     0/1:43:1  
      ...

    $ cat mutect2.03.merge.sitesOnly.vcf
      #CHROM  POS ID  REF  ALT  QUAL FILTER  INFO               
      chrM    64  .   C    T    .    .       AC=1;AN=2;HV       
      chrM    73  .   A    G    .    .       AC=1;AN=2;HV;NUMT  
      chrM    73  .   A    G    .    .       AC=2;AN=4;HV       
      ....
     
#### SNV counts ####

    $ cat mutect2.03.tab            
      Run                        H           h   S   s   I  i  Hp  hp  Sp  sp  Ip  ip  
      chrM.A                     31          44  28  36  3  8  3   0   0   0   3   0   
      chrM.B                     27          42  25  34  2  8  5   0   3   0   2   0   
      chrM.C                     39          43  35  35  4  8  2   0   0   0   2   0   

#### Haplogroups ####

    $ cat mutect2.haplogroup.tab    
      Run                        haplogroup  
      chrM.A                     A2+(64)     
      chrM.B                     B2          
      chrM.C                     C           

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

## EXAMPLE 3 : CUSTOM FILTERING ##

### using init.sh ###

    $ nano init.sh
      export FNAME="no_qual_haplotype_strand"
      export FRULE="egrep -v \"qual|haplotype|strand\""
    # rerun run.sh ...

### command line ###    

    $ M=mutect2
    $ cat $ODIR/$M.00.concat.vcf  | egrep -v 'qual|haplotype|strand'  > \
        $ODIR/$M.no_qual_haplotype_strand.00.concat.vcf
    $ $SDIR/snpCount.sh $IN $ODIR $M.no_qual_haplotype_strand 03
    $ $SDIR/snpCount.sh $IN $ODIR $M.no_qual_haplotype_strand 05
    $ $SDIR/snpCount.sh $IN $ODIR $M.no_qual_haplotype_strand 10

