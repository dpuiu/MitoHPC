# HP : Mitochondrial Heteroplasmy Pipeline # 

## INPUT ##
 
    * Illumina paired-end reads pre-aligned to a reference: .bam or .cram files 
    * Alignment indeces: .bai or .crai files
    * Alignment indexstats: .idxstats files

## SYSTEM PREREQUISITES ##

    * OS:                UNIX/LINUX 
    * SOFTWARE PACKAGES: git,wget,tar,unzip,make,autoconf,gcc,java,perl,python 
 
## PIPELINE PREREQUISITES ##

    * SOFTWARE PACKAGES: bwa,samtools,bedtools,fastp,samblaster,bcftools,htslib,freebayes
    * JAVA JARS:         gatk,mutserve,haplogrep,haplocheck
    * HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    $ git clone https://github.com/dpuiu/HP.git

### UPDATE PIPELINE (optional) ###

    $ cd HP/
    $ git pull
      or
    $ git checkout .

### SETUP ENVIRONMENT ###

    $ cd HP/scripts
    $ export HP_SDIR=`pwd`                                       # set script directory variable 
    $ . ./init.sh                                                # or init.hs38DH.sh or init.mm39.sh

### INSTALL SYSTEM PREREQUISITES (optional) ###

    $ sudo $HP_SDIR/install_sysprerequisites.sh                  # install perl,pthon,java,wget ...

### INSTALL PREREQUISITES ; CHECK INSTALL ###

    $ $HP_SDIR/install_prerequisites.sh                          # bwa,samtools,bedtools ...
      or
    $ $HP_SDIR/install_prerequisites.sh -f                       # force install

    $ $HP_SDIR/checkInstall.sh
    # if successfull => "Success message!"

    $ cat checkInstall.log                                       # records software version/path info

## PIPELINE USAGE ##

### SETUP ENVIRONMENT ###

    # go to your working directory ; copy over the init.sh file
    $ cp -i $HP_SDIR/init.sh .                                   # copy init to working dir.


    # edit init.sh file
    $ nano init.sh                                               # check/edit local init file
        ...
       export HP_ADIR=$PWD/bams/                                 # alignment dir; .bam, .bai, [.idxstats] files
       # or  
       export HP_ADIR=$PWD/crams/                                # alignment dir; .cram, .crai, [.idxstats] files
      
       export HP_ODIR=$PWD/out/                                  # output dir  
       export HP_IN=$PWD/in.txt                                  # input file

       export HP_L=222000                                        #  Use at most 222K reads

       export HP_DOPT="--removeDups"                             # samblaster deduplication option
       export HP_GOPT=                                           # GATK mutect2 options; Ex: "-max-reads-per-alignment-start 20 -mitochondria-mode"  
       export HP_FOPT=                                           # FASTP read trimming options; Ex: "-q 20 -e 20"              

       export HP_FRULE=                                          # VCF output filtering options; Ex: "bcftools view -f PASS,clustered_events,multiallelic"

       find $HP_ADIR/ -name "*.bam" -o -name "*.cram" | \
        $HP_SDIR/ls2in.pl -out $HP_ODIR | sort -V > $HP_IN       # generate input file; 3 column tsv file
								 # prefix, input file(bam/cram), outpt prefix
      ...                                                        # can be edited

       export HP_SH=bash                                         # job scheduling: bash, qsub,sbatch, ..

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

    # under $HP_ODIR: TAB/SUMMARY/VCF/FASTA Files: 

    count.tab                                                    # total reads & mtDNA-CN counts
    subsample.tab                                                # subsampled read counts
    cvg.tab                                                      # subsampled coverage stats

    # 1st ITERATION
    {mutect2,mutserve,freebayes}.{03,05,10}.{concat,merge[.sitesOnly]}.vcf # SNVs; 3,5,10% heteroplasmy thold
    {mutect2,mutserve,freebayes}.{03,05,10}.tab                            # SNV counts
    {mutect2,mutserve,freebayes}.{03,05,10}.summary                        # SNV count summaries

    {mutect2,mutserve,freebayes}.fa                                        # consensus sequence
    {mutect2,mutserve,freebayes}.haplogroup[1].tab                         # haplogroup
    {mutect2,mutserve,freebayes}.haplocheck.tab                            # contamination screen   

    # 2nd ITERATION
    {mutect2_mutect2.mutect2,freebayes_freebayes.freebayes}.{03,05,10}.{concat,merge[.sitesOnly]}.vcf	
    {mutect2_mutect2.mutect2,freebayes_freebayes.freebayes}.{03,05,10}.{tab,summary}

## EXAMPLE 1 ##

    3 simulated samples, 100x chrM coverage 

### INPUT ###

    $ head $HP_IN
      #sampleName     inputFile          outputPath/prefix
      chrM.A          bams/chrM.A.bam    out/chrM.A/chrM.A
      chrM.B          bams/chrM.B.bam    out/chrM.B/chrM.B
      chrM.C          bams/chrM.C.bam    out/chrM.C/chrM.C
      ...

### OUTPUT ###

#### output directories : 1/sample ####

    $ ls out/
      chrM.A/
      chrM.B/
      chrM.C/
      ...

#### directory structure ####

    $ ls $HP_ADIR/*
      bams/chrM.A.bam
      bams/chrM.A.bam.bai
      bams/chrM.A.idxstats
      bams/chrM.A.count
     ...

#### read counts and mtDNA-CN(M) ####

     $ head $HP_ODIR/count.tab
       Run       all        mapped     chrM    M
       chrM.A    851537886  848029490  396766  182
       chrM.B    884383716  882213718  506597  223
       chrM.C    786560467  785208588  503241  249
       ...

#### vcf files : concat & merge ####

    # 1st iteration : homoplasmies(AF=1) and heteroplasmies(AF<1)
    $ cat mutect2.03.concat.vcf  
      ...
      ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="INDEL">
      ##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">
      ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
      ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
      ##INFO=<ID=HP,Number=0,Type=Flag,Description="Homoloplymer">
      ##INFO=<ID=HS,Number=0,Type=Flag,Description="Hot spot">
      ##INFO=<ID=CDS,Number=0,Type=Flag,Description="coding region">
      ##INFO=<ID=RNR,Number=0,Type=Flag,Description="rRNA">
      ##INFO=<ID=TRN,Number=0,Type=Flag,Description="tRNA">
      ##INFO=<ID=DLOOP,Number=0,Type=Flag,Description="D-loop">
      ##INFO=<ID=HG,Number=1,Type=String,Description="haplogroup">
      ##INFO=<ID=NUMT,Number=0,Type=Flag,Description="NUMT like">
      ##INFO=<ID=HV,Number=0,Type=Flag,Description="Hypervariable">
      ##INFO=<ID=AP,Number=1,Type=String,Description="Mitimpact APOGEE">
      ##INFO=<ID=APS,Number=1,Type=String,Description="Mitimpact APOGEE_score">
      ...
      #CHROM POS   ID REF ALT  QUAL  FILTER                      INFO                                                                               FORMAT SAMPLE
      chrM   64    .  C   T    .     clustered_events;haplotype  HV=HV-II;DLOOP;GT=0|1;DP=66;AF=1                                                   SM     chrM.A
      chrM   73    .  A   G    .     PASS                        HV=HV-II;DLOOP;NUMT=MT780574.1|KM281532.1|KM281533.1|KM281524.1;GT=0/1;DP=51;AF=1  SM     chrM.B
      chrM   73    .  A   G    .     PASS                        HV=HV-II;DLOOP;NUMT=MT780574.1|KM281532.1|KM281533.1|KM281524.1;GT=0/1;DP=46;AF=1  SM     chrM.C
      chrM   73    .  A   G    .     clustered_events;haplotype  HV=HV-II;DLOOP;NUMT=MT780574.1|KM281532.1|KM281533.1|KM281524.1;GT=0|1;DP=70;AF=1  SM     chrM.A
      chrM   146   .  T   C    .     clustered_events;haplotype  HV=HV-II;DLOOP;GT=0|1;DP=88;AF=1                                                   SM     chrM.A
      ...
      chrM   1037  .  A   C    .     PASS                        RNR=RNR1;GT=0/1;DP=70;AF=0.178                                                     SM     chrM.B
      chrM   1038  .  C   G    .     PASS                        RNR=RNR1;GT=0/1;DP=70;AF=0.165                                                     SM     chrM.A
      chrM   1042  .  T   A    .     PASS                        RNR=RNR1;GT=0/1;DP=60;AF=0.24                                                      SM     chrM.C
      ...
      chrM   3988  .  T   G    .     PASS                        CDS=ND1;AP=Pathogenic;APS=0.52;GT=0/1;DP=60;AF=0.224;AP=Pathogenic;APS=0.52        SM     chrM.B 
      chrM   3989  .  A   T    .     PASS                        CDS=ND1;AP=Pathogenic;APS=0.52;GT=0/1;DP=65;AF=0.237;AP=Pathogenic;APS=0.52        SM     chrM.A
      chrM   5186  .  A   T    .     PASS                        CDS=ND2;HG=U;AP=Pathogenic;APS=0.51;GT=0|1;DP=97;AF=0.212;AP=Pathogenic;APS=0.51   SM     chrM.C
      ...

    $ cat mutect2.03.merge.vcf  
      ...      
      #CHROM POS   ID REF ALT  QUAL  FILTER  INFO                                      FORMAT    chrM.A        chrM.B        chrM.C
      ...
      chrM   3988  .  T   G    .     .       AC=1;AN=2;CDS=ND1;AP=Pathogenic;APS=0.52  GT:DP:AF  .:.:.         0/1:60:0.224  .:.:.
      chrM   3989  .  A   T    .     .       AC=1;AN=2;CDS=ND1;AP=Pathogenic;APS=0.52  GT:DP:AF  0/1:65:0.237  .:.:.         .:.:.
      chrM   3993  .  A   T    .     .       AC=1;AN=2;CDS=ND1                         GT:DP:AF  .:.:.         .:.:.         0/1:52:0.258
      ...

    $ cat mutect2.03.merge.sitesOnly.vcf
       ...
      #CHROM POS   ID REF ALT  QUAL  FILTER  INFO                                      
      ...
      chrM   3988  .  T   G    .     .       AC=1;AN=2;CDS=ND1;AP=Pathogenic;APS=0.52
      chrM   3989  .  A   T    .     .       AC=1;AN=2;CDS=ND1;AP=Pathogenic;APS=0.52
      chrM   3993  .  A   T    .     .       AC=1;AN=2;CDS=ND1                       
      ...
     
#### SNV counts ####

    # 1st iteration
    $ cat mutect2.03.tab            
      Run     H   h   S   s   I  i  Hp  hp  Sp  sp  Ip ip  A
      chrM.A  30  45  28  36  2  9  28  44  28  36  0  8   75
      chrM.B  26  43  25  34  1  9  22  42  22  34  0  8   69
      chrM.C  39  43  35  35  4  8  37  43  35  35  2  8   82
      ...

    # 2nd iteration
    $ cat  mutect2.mutect2.03.tab 
      Run     H   h   S   s   I  i  Hp  hp  Sp  sp  Ip  ip  A
      chrM.A  30  44  28  36  2  8  28  44  28  36  0   8   74		# one less heteroplasmy (h) compared with the 1st iteration
      chrM.B  26  42  25  34  1  8  22  42  22  34  0   8   68		# one less heteroplasmy	(h) compared with the 1st iteration
      chrM.C  39  43  35  35  4  8  37  43  35  35  2   8   82
      ...

#### SNV Summaries ####

    # 1st iteration
    $ cat mutect2.03.summary  | column -t
      id  count  nonZero  min  max  median  mean   sum
      H   3      3        26   39   30      31.67  95
      h   3      3        43   45   43      43.67  131
      S   3      3        25   35   28      29.33  88
      s   3      3        34   36   35      35     105
      I   3      3        1    4    2       2.33   7
      i   3      3        8    9    9       8.67   26
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

#### Haplocheck Contamination ####

    $ cat mutect2.haplocheck.tab  | column -t
      Run     ContaminationStatus  ContaminationLevel  Distance  SampleCoverage
      chrM.A  NO                   ND                  14        78
      chrM.B  NO                   ND                  14        79
      chrM.C  NO                   ND                  14        78
      ...

#### Custom Filtering ####

    # define filter name and command
    $ nano init.sh
      export HP_FNAME="no_qual_haplotype_strand"
      export HP_FRULE="egrep -v \"qual|haplotype|strand\""

    # reinit and rerun
    $ . ./init.sh 
    $ $HP_SDIR/run.sh > run.all.sh 
    $ bash ./run.all.sh  

    # check additional output files
    $ ls $HP_OUT/$HP_M.no_qual_haplotype_strand.*

### LEGEND ###

    Run                           # Sample name

    all                           # number of reads in the sample reads
    mapped                        # number of aligned reads in the sample
    MT                            # number of reads aligned to chrM
    M                             # mtDNA copy number

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
