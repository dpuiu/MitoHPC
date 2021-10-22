# HP : Heteroplasmy Pipeline # 


## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINIX 
    SOFTWARE PACKAGES: git, wget, tar, unzip, make, autoconf, gcc(v4.8.5+), java(v1.8.0+), perl, python 
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa, samtools, bedtools, fastp, samblaster, bcftools, htslib, vcftools, st
    JAVA JARS:         gatk, mutserve, haplogrep, haplocheck
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/HP.git

### SETUP ENVIRONMENT (1) ###
    
    $ PWD=`pwd`
    $ export HP_SDIR=path_to_pipeline_home_directory/HP/scripts/   # set script directory variable

### INSTALL PIPELINE PREREQUISITES ; CHECK INSTALL ###

    $ sudo $HP_SDIR/install_sysprerequisites.sh
    $ $HP_SDIR/install_prerequisites.sh  
      
    # if successfull => "Success message!"

    $ $HP_SDIR/checkInstall.sh
    $ cat checkInstall.log

## USAGE ##

### COPY init.sh SCRIPT ###

    # go to the working directory ; copy init.sh 

    $ cp -i $HP_SDIR/init.sh .

### SETUP ENVIRONMENT (2) ###

    $ nano init.sh 
        ...
       export HP_ADIR=$PWD/bams/                                           # alignment dir
       export HP_ODIR=$PWD/out/                                            # output dir  
       export HP_IN=$PWD/in.txt                                            # input file
    
    $ . ./init.sh                                                          # source init file 

    $ printenv | grep HP_                                                  # check variables 

    $ mkdir -p $HP_ODIR                                       

    $ find $HP_ADIR/ -name "*.bam" -o -name "*.cram" | \
        $HP_SDIR/ls2in.pl -out $HP_ODIR | sort > $HP_IN                    # generate input file

### GENERATE ALIGNMENT INDEX AND READ COUNT FILES; COMPUTE mtDNA-CN ####

     $ cut -f2 $HP_IN  | sed "s|^|$HP_SH $HP_SDIR/samtools.sh |" > samtools.all.sh

     $ bash ./samtools.all.sh

     # wait for the jobs to finish ...

     $ cut -f2 $HP_IN | sed -r 's|(.*)\.|\1\t|g' | cut -f1 | \
         sed 's|$|.count|' | xargs cat | \
         $HP_SDIR/uniq.pl | $HP_SDIR/getCN.pl > $HP_ODIR/count.tab

### SUBSAMPLE ALIGNMENT FILES (optional) ###

        $ nano init.sh
          ...
          export HP_L=222000                                                # number of reads
          export HP_LDIR=...                                                # subsample dir

       $ join $HP_IN $HP_ODIR/count.tab | \                                 
           perl -ane '$s=$ENV{HP_L}/$F[5]; $s=($s<1)?"-s $s":""; \
             print "samtools view $s $F[1] $ENV{HP_MT} $ENV{HP_NUMT} -T $ENV{HP_RDIR}/$ENV{HP_H}.fa -b \
             $ENV{HP_LDIR}/$F[0].bam; \
             samtools index $ENV{HP_LDIR}/$F[0].bam\n"' | bash              # subsample

       $ find $HP_LDIR/ -name "*.bam" -o -name "*.cram" | \                 
          $HP_SDIR/ls2in.pl -out $HP_ODIR | sort > $HP_IN		    # re-create input file
  
### RUN PIPELINE  ###
 
    $ $HP_SDIR/run.sh > filter.all.sh

    $ bash ./filter.all.sh                                             


### RE-RUN PIPELINE (optional) ###

    $ nano init.sh           						    # update parameters if needed                                               

    $ $HP_SDIR/run.sh > filter.all.sh                                       

    $ bash ./filter.all.sh
     
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

    $ head $HP_IN
      #sampleName     inputFile          outputPath/prefix
      chrM.A          bams/chrM.A.bam    out/chrM.A/chrM.A
      chrM.B          bams/chrM.B.bam    out/chrM.B/chrM.B
      chrM.C          bams/chrM.C.bam    out/chrM.C/chrM.C
      ...

 ### read counts and mtDNA-CN(M)

     $ head $HP_ODIR/count.tab
       Run	 all        mapped     chrM    M
       chrM.A    851537886  848029490  396766  181.7
       chrM.B    884383716  882213718  506597  223.01
       chrM.C    786560467  785208588  503241  248.9
       ...


### after running samtools.sh ###
    $ ls $HP_ADIR/*
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

##### 1st itteration #####

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
      Run     H   h   S   s   I  i  Hp  hp  Sp  sp  Ip  ip  A
      chrM.A  31  43  28  35  3  8  28  43  28  35  0   8   74
      chrM.B  27  43  25  35  2  8  22  41  22  33  0   8   70
      chrM.C  40  42  36  34  4  8  38  41  36  33  2   8   82
      ...

    $ cat  mutect2.mutect2.03.tab 
      Run     H  h   S  s   I  i  Hp  hp  Sp  sp  Ip  ip  A
      chrM.A  0  43  0  35  0  8  0   43  0   35  0   8   43
      chrM.B  0  43  0  35  0  8  0   41  0   33  0   8   43
      chrM.C  0  42  0  34  0  8  0   41  0   33  0   8   42

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

    H                             # number of homozygous SNVs (the multiallelic SNVs counted once)
    h                             # number of heterozygous SNVs
    S                             # number of homozygous SNVs (excluding INDELs)
    s                             # number of heterozygous SNVs (excluding INDELs)
    I                             # number of homozygous INDELs
    i                             # number of heterozygous INDELs
    ...
    ?p                            # "p" suffix stands for non homopolimeric
    ...
    A                             # total number of SNVs (H+h)

## EXAMPLE 3 : CUSTOM FILTERING ##

### using init.sh ###

    $ nano init.sh
      export HP_FNAME="no_qual_haplotype_strand"
      export HP_FRULE="egrep -v \"qual|haplotype|strand\""
    # rerun run.sh ...

### command line ###    

    $ HP_M=mutect2
    $ cat $HP_ODIR/$M.00.concat.vcf  | egrep -v 'qual|haplotype|strand'  > $HP_ODIR/$M.no_qual_haplotype_strand.00.concat.vcf
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 03
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 05
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 10

