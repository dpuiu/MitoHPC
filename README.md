# HP : Mitochondrial Heteroplasmy Pipeline # 

## Input: Illumina paired-end reads pre-aligned to a reference ##

## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINUX 
    SOFTWARE PACKAGES: git,wget,tar,unzip,make,autoconf,gcc,java,perl,python 
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa,samtools,bedtools,fastp,samblaster,bcftools,htslib
    JAVA JARS:         gatk,mutserve,haplogrep,haplocheck
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/HP.git

### UPDATE PIPELINE (optional) ###

    $ cd HP/
    $ git pull

### SETUP ENVIRONMENT ###

    $ cd HP/scripts
    $ export HP_SDIR=`pwd`                                       # set script directory variable 

### INSTALL SYSTEM PREREQUISITES (optional) ###

    $ sudo $HP_SDIR/install_sysprerequisites.sh                  # install perl,pthon,java,wget ...

### INSTALL PREREQUISITES ; CHECK INSTALL ###

    $ $HP_SDIR/install_prerequisites.sh                          # bwa,samtools,bedtools ...
      or
    $ $HP_SDIR/install_prerequisites.sh -f                       # force install

    $ $HP_SDIR/checkInstall.sh
    # if successfull => "Success message!"

    $ cat checkInstall.log                                       # records software version/path info

## USAGE ##

### SETUP ENVIRONMENT ###

    # go to your working directory ; copy over the init.sh file
    $ cp -i $HP_SDIR/init.sh .                                   # copy init to working dir.

    # edit init.sh file
    $ nano init.sh                                               # check/edit local init file
        ...
       export HP_ADIR=$PWD/bams/                                 # alignment dir; contains .bam, .bam.bai and/or .idxstats file
       # or  
       export HP_ADIR=$PWD/crams/                                 # alignment dir; contains .cram, .cram.crai and/or .idxstats file
      
       export HP_ODIR=$PWD/out/                                  # output dir  
       export HP_IN=$PWD/in.txt                                  # input file

       export HP_L=222000                                        #  Use at most 222K reads

       export HP_DOPT="--removeDups"                             # samblaster deduplication option
       export HP_GOPT="-max-reads-per-alignment-start 20 -mitochondria-mode"  # GATK mutect2 options
       export HP_FOPT="-q 20 -e 20"                             # FASTP options

       find $HP_ADIR/ -name "*.bam" -o -name "*.cram" | \
        $HP_SDIR/ls2in.pl -out $HP_ODIR | sort -V > $HP_IN       # generate input file; 3 column
      ...                                                        # tab delimited, can be edited

    $ . ./init.sh                                                # source init file 

    $ printenv | grep HP_ | sort                                 # optional ; check HP_ variables 
    $ nano $HP_IN                                                # optional ; check input file

### RUN PIPELINE  ###
 
    $ $HP_SDIR/run.sh > run.all.sh                               # create command file in working dir.
    $ bash ./run.all.sh                                          # execute command file      

### RE-RUN PIPELINE (optional) ###

    $ nano init.sh                                               # update parameters
    $ . ./init.sh                                                # source init file
    $ $HP_SDIR/run.sh > run.all.sh                               # recreate command file

    $ bash ./run.all.sh
     
## OUTPUT ##

    # under $HP_ODIR:

    TAB/SUMMARY/VCF/FASTA Files: 

    all.count.tab                                                # total reads & mtDNA-CN counts
    count.tab                                                    # subsampled read counts
    cvg.tab                                                      # subsampled coverage stats

### 1st ITERATION ### 

    {mutect2,mutserve}.{03,05,10}.{concat,merge[.sitesOnly]}.vcf # SNVs; 3,5,10% heteroplasmy thold
    {mutect2,mutserve}.{03,05,10}.tab                            # SNV counts
    {mutect2,mutserve}.{03,05,10}.summary                        # SNV count summaries

    {mutect2,mutserve}.fa                                        # consensus sequence
    {mutect2,mutserve}.haplogroup[1].tab                         # haplogroup
    {mutect2,mutserve}.haplocheck.tab                            # contamination screen   


### 2nd ITERATION ###

    mutect2.mutect2.{03,05,10}.{concat,merge[.sitesOnly]}.vcf	
    mutect2.mutect2.{03,05,10}.tab

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

### directory structure ###

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

    # 1st iteration : homoplasmies(AF=1) and heteroplasmies(AF<1)
    $ cat mutect2.03.concat.vcf  
      ...
      ##INFO=<ID=SM,Number=1,Type=String,Description="Sample">
      ##INFO=<ID=HP,Number=0,Type=Flag,Description="Homoloplymer">
      ##INFO=<ID=HS,Number=0,Type=Flag,Description="Hot spot">
      ##INFO=<ID=CDS,Number=0,Type=Flag,Description="coding region">
      ##INFO=<ID=RNR,Number=0,Type=Flag,Description="rRNA">
      ##INFO=<ID=TRN,Number=0,Type=Flag,Description="tRNA">
      ##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="D-loop">
      ##INFO=<ID=HG,Number=1,Type=String,Description="haplogroup">
      ##INFO=<ID=NUMT,Number=0,Type=Flag,Description="NUMT like">
      ##INFO=<ID=HV,Number=0,Type=Flag,Description="Hypervariable">
      ##INFO=<ID=SM,Number=1,Type=String,Description="Sample">
      ##INFO=<ID=CR,Number=1,Type=String,Description="CADD Raw score">
      ##INFO=<ID=CP,Number=1,Type=String,Description="CADD PHREDD score">
      ...
      ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
      ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
      ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
      ...
      #CHROM  POS  ID  REF  ALT  QUAL  FILTER                        INFO                                          FORMAT    SAMPLE
      chrM    64   .   C    T    .     clustered_events;haplotype    SM=chrM.A;HV;DLOOP;CR=1.140555;CP=11.43       GT:DP:AF  0|1:67:1
      chrM    73   .   A    G    .     PASS                          SM=chrM.B;HV;DLOOP;NUMT;CR=1.039800;CP=10.89  GT:DP:AF  0/1:52:1
      chrM    73   .   A    G    .     PASS                          SM=chrM.C;HV;DLOOP;NUMT;CR=1.039800;CP=10.89  GT:DP:AF  0/1:43:1
      chrM    73   .   A    G    .     clustered_events;haplotype    SM=chrM.A;HV;DLOOP;NUMT;CR=1.039800;CP=10.89  GT:DP:AF  0|1:70:1
      chrM    146  .   T    C    .     clustered_events;haplotype    SM=chrM.A;HV;DLOOP;CR=0.949142;CP=10.37       GT:DP:AF  0|1:88:1
      ...
      chrM    374  .   A    G    .     clustered_events              SM=chrM.A;DLOOP;CR=0.998787;CP=10.66          GT:DP:AF  0/1:76:0.139
      chrM    375  .   C    T    .     clustered_events;strand_bias  SM=chrM.B;DLOOP;CR=1.074472;CP=11.08          GT:DP:AF  0/1:65:0.162
      chrM    378  .     C    T    .   clustered_events;strand_bias  SM=chrM.C;DLOOP;CR=1.043995;CP=10.91          GT:DP:AF  0/1:82:0.128
      ...
      chrM    1037 .   A    C    .     PASS                          SM=chrM.B;RNR;CR=0.682748;CP=8.714            GT:DP:AF  0/1:70:0.178
      chrM    1038 .   C    G    .     PASS                          SM=chrM.A;RNR;CR=0.822718;CP=9.603            GT:DP:AF  0/1:70:0.165
      chrM    1042 .   T    A     .    PASS                          SM=chrM.C;RNR;CR=1.104043;CP=11.24            GT:DP:AF  0/1:60:0.24
      ...
    $ cat mutect2.03.merge.vcf  
      ...  
      #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO                                          FORMAT    chrM.A    chrM.B    chrM.C
      chrM    64   .   C    T    .     .       AC=1;AN=2;HV;DLOOP;CR=1.140555;CP=11.43       GT:DP:AF  0|1:67:1  .:.:.     .:.:.
      chrM    73   .   A    G    .     .       AC=3;AN=6;HV;DLOOP;NUMT;CR=1.039800;CP=10.89  GT:DP:AF  0|1:70:1  0/1:52:1  0/1:43:1
      chrM    146  .   T    C    .     .       AC=1;AN=2;HV;DLOOP;CR=0.949142;CP=10.37       GT:DP:AF  0|1:88:1  .:.:.     .:.:.
      ...

    $ cat mutect2.03.merge.sitesOnly.vcf
     ...
      #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
      chrM    64   .   C    T    .     .       AC=1;AN=2;HV;DLOOP;CR=1.140555;CP=11.43
      chrM    73   .   A    G    .     .       AC=3;AN=6;HV;DLOOP;NUMT;CR=1.039800;CP=10.89
      chrM    146  .   T    C    .     .       AC=1;AN=2;HV;DLOOP;CR=0.949142;CP=10.37
      ...

    # 2nd iteration ; mostly heteroplasmies(AF<1)
    $ cat mutect2.mutect2.03.concat.vcf 
      #CHROM  POS   ID  REF   ALT  QUAL  FILTER  INFO                                  FORMAT    SAMPLE
      chrM    374   .   A     G    .     PASS    SM=chrM.A;DLOOP;CR=0.998787;CP=10.66  GT:DP:AF  0/1:62:0.138
      chrM    375   .   C     T    .     PASS    SM=chrM.B;DLOOP;CR=1.074472;CP=11.08  GT:DP:AF  0/1:48:0.178
      chrM    378   .   C     T    .     PASS    SM=chrM.C;DLOOP;CR=1.043995;CP=10.91  GT:DP:AF  0/1:71:0.121
      ...
      chrM    1037  .   A     C    .     PASS    SM=chrM.B;RNR;CR=0.682748;CP=8.714    GT:DP:AF  0/1:70:0.178
      chrM    1038  .   C     G    .     PASS    SM=chrM.A;RNR;CR=0.822718;CP=9.603    GT:DP:AF  0/1:70:0.165
      chrM    1042  .   T     A    .     PASS    SM=chrM.C;RNR;CR=1.104043;CP=11.24    GT:DP:AF  0/1:60:0.24
      ...
     
#### SNV counts ####

    # 1st iteration
    $ cat mutect2.03.tab            
      Run     H   h   S   s   I  i  Hp  hp  Sp  sp  Ip  ip  A
      chrM.A  31  43  28  35  3  8  28  43  28  35  0   8   74
      chrM.B  27  43  25  35  2  8  22  41  22  33  0   8   70
      chrM.C  40  42  36  34  4  8  38  41  36  33  2   8   82
      ...

    # 2nd iteration ; no homoplamies(H=S=I=0)
    $ cat  mutect2.mutect2.03.tab 
      Run     H  h   S  s   I  i  Hp  hp  Sp  sp  Ip  ip  A
      chrM.A  0  43  0  35  0  8  0   43  0   35  0   8   43
      chrM.B  0  43  0  35  0  8  0   41  0   33  0   8   43
      chrM.C  0  42  0  34  0  8  0   41  0   33  0   8   42
      ...

#### Summaries ####
    # 1st iteration
    $ cat mutect2.03.summary  | column -t
      id  count  nonZero  min  max  median  mean   sum
       H   3      3        27   39   31      32.33  97
       h   3      3        42   44   43      43     129
       S   3      3        25   35   28      29.33  88
       s   3      3        34   36   35      35     105
       I   3      3        2    4    3       3      9
       i   3      3        8    8    8       8      24
       Hp  3      3        22   37   28      29     87
       hp  3      3        42   44   43      43     129
       Sp  3      3        22   35   28      28.33  85
       sp  3      3        34   36   35      35     105
       Ip  3      1        0    2    0       0.67   2
       ip  3      3        8    8    8       8      24
       A   3      3        69   82   75      75.33  226

#### Haplogroups ####

    $ cat mutect2.haplogroup.tab    
      Run                        haplogroup  
      chrM.A                     A2+(64)     
      chrM.B                     B2          
      chrM.C                     C           
      ...

#### Contamination ####

    $ cat mutect2.haplocheck.tab  | column -t
      Run     ContaminationStatus  ContaminationLevel  Distance  SampleCoverage
      chrM.A  NO                   ND                  14        78
      chrM.B  NO                   ND                  14        79
      chrM.C  NO                   ND                  14        78
      ...

## EXAMPLE 3 : CUSTOM FILTERING ##

### using init.sh ###

    $ nano init.sh
      export HP_FNAME="no_qual_haplotype_strand"
      export HP_FRULE="egrep -v \"qual|haplotype|strand\""
    # rerun run.sh ...

### command line ###    

    $ HP_M=mutect2
    $ cat $HP_ODIR/$M.00.concat.vcf  | egrep -v 'qual|haplotype|strand'  > \
        $HP_ODIR/$M.no_qual_haplotype_strand.00.concat.vcf
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 03
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 05
    $ $HP_SDIR/snpCount.sh $HP_IN $ODIR $HP_M.no_qual_haplotype_strand 10

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
