# HP : Heteroplasmy Pipeline # 

## SYSTEM PREREQUISITES ##

    OS:                UNIX/LINIX 
    SOFTWARE PACKAGES: git, wget, tar, unzip, autoconf, gcc(v4.8.5+), java(v1.8.0+), perl(v5.16.3+) ....
 
## PIPELINE PREREQUISITES ##

    SOFTWARE PACKAGES: bwa, samtools, bedtools, fastp, samblaster, bcftools, htslib, vcftools, st
    JAVA JARS:         picard, mutserve, gatk, haplogrep
    HUMAN ASSEMBLY:    hs38DH

## INSTALL ## 

### DOWNLOAD PIPELINE ###

    git clone https://github.com/dpuiu/HP.git

### INSTALL SYSTEM PREREQUISITES (optional) ###

    $SDIR/install_sysprerequisites.sh

### SETUP ENVIRONMENT ###
    
    export SDIR=$PWD/HP/scripts/   # set SDIR variable (script directory)
    echo $SDIR                     # check SDIR variable is set correctly
 
    source $SDIR/init.sh
    ... or    
    source $SDIR/init_marcc.sh     # for MARCC computer cluster, uses SLURM job scheduler

    # could be added to ~/.bashrc             # set environmental variables
    echo "## HP settings ##"  >> ~/.bashrc
    echo "SDIR=$SDIR" >> ~/.bashrc
    echo "source $SDIR/init.sh" >> ~/.bashrc

### INSTALL PIPELINE PREREQUISITES (optional) ###

    $SDIR/install_prerequisites.sh  # install PIPELINE PREREQUISITES

### CHECK INSTALL ###
  
    # if successfull => "Success message!"
    $SDIR/checkInstall.sh
    cat checkInstall.log

## USAGE ##

### GENERATE INPUT FILE  ###

    find bams/ -name "*.bam"  > in.txt
    ... or
    find crams/-name "*.cram" > in.txt

### SPLIT INPUT FILE (optional; for large datasets) ###
   
    # Ex: up to 9 sets of 100

    split -d -a 1 --numeric=1 -l 100 in.txt  in. --additional-suffix=.txt
    ... or
    split -d -a 2 --numeric=1 -l 100 in.txt  in. --additional-suffix=.txt
   
### GENERATE INDEX AND COUNT FILES ###

     # prefix.bam  => prefix.bam.bai & prefix.count 
     .... or 
     # prefix.cram => prefix.cram.crai & prefix.count 

     sed 's|^|samtools.sh |' in.txt > samtools.all.sh
     ... or (MARCC)
     sed 's|^|sbatch --p shared --time=24:0:0 samtools.sh|' in.txt > samtools.all.sh
     chmod a+x ./samtools.all.sh
     ./samtools.all.sh

### GENERATE PIPLEINE SCRIPT ###

    $SDIR/run.sh in.txt out/  > filter.all.sh
    ... or
    $SDIR/run.sh in.1.txt out/ > filter.1.sh
    $SDIR/run.sh in.2.txt out/ > filter.2.sh
    $SDIR/run.sh in.3.txt out/ > filter.3.sh
    ...
    chmod u+x ./filter.*sh

### RUN PIPELINE  ###

    nohup ./filter.all.sh &
    ... or
    nohup ./filter.1.sh &
    nohup ./filter.2.sh &
    nohup ./filter.3.sh &
    ...
    ls filter.?.sh | parallel
    ... or (MARCC)
    sbatch --p shared --time=24:0:0 ./filter.all.sh

## OUTPUT ##

### TAB/VCF Files: ###

    # 1st itteration (ref=rCRS or RSRS, snp caller mutect2 or mutserve)
    count.tab 
 
    mutect2.03.vcf	# HOM & HET SNPs; 3% min HET THOLD
    mutect2.05.vcf	#                 5%
    mutect2.10.vcf	#                 10%
    count_mutect2.tab   # read, MT & SNP counts ; SNPs: HOM & HET 3,5,10% thold 
    ... or
    mutserve.03.vcf
    mutserve.05.vcf
    mutserve.10.vcf
    count_mutserve.tab

    # optionl ; if mutect2 is run a 2nd time, ref=sample MT consensus
    mutect2.mutect2.03.vcf
    mutect2.mutect2.05.vcf
    mutect2.mutect2.10.vcf
    count_mutect2.mutect2.tab

## EXAMPLES ##

### use RSRS for realignment ###

    #script      inFile outDir   humanRef  MTRef
    $SDIR/run.sh in.txt out/    
    ... or
    $SDIR/run.sh in.txt out/     hs38DH.fa RSRS.fa

### use rCRS for realignment, mutserve for SNP calling ###

    $SDIR/run.sh in.txt out/ hs38DH.fa rCRS.fa mutserve

### output ; simulated haplogroup A,B,C datasets ###

    head count.tab 
      Run     all        mapped     chrM    filter  M
      chrM.A  740589366  739237125  487382  205095  256.05
      chrM.B  763658318  762297733  495743  205156  252.56
      chrM.C  749938586  748667901  590963  200121  306.55

    head mutect2.03.vcf 
      ##fileformat=VCFv4.2
      ##source=Mutect2
      ##reference=file://../..//RefSeq//rCRS.fa>
      ##contig=<ID=rCRS,length=16569>
      #CHROM	POS	ID	REF	ALT	QUAL	FILTER	                        INFO	                FORMAT	SAMPLE
      rCRS	64	.	C	T	.	clustered_events;haplotype	SNP;DP=1290	        SM	chrM.A
      rCRS	73	.	A	G	.	clustered_events;haplotype	SNP;DP=1421	        SM	chrM.A
      rCRS	73	.	A	G	.	PASS	                        SNP;DP=1002	        SM	chrM.B
      rCRS	73	.	A	G	.	PASS	                        SNP;DP=833	        SM	chrM.C
      rCRS	146	.	T	C	.	clustered_events;haplotype	SNP;DP=1684	        SM	chrM.A
      rCRS	153	.	A	G	.	clustered_events;haplotype	SNP;DP=1695	        SM	chrM.A
      ...
      rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=1416;HP	SM	chrM.A
      rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=1227;HP	SM	chrM.B
      rCRS	310	.	T	TC	.	clustered_events	        INDEL;DP=1978;HP	SM	chrM.C
      ...
      rCRS	374	.	A	G	.	clustered_events	        SNP;DP=1469;AF=0.139	SM	chrM.A
      rCRS	375	.	C	T	.	clustered_events;strand_bias	SNP;DP=1250;AF=0.162	SM	chrM.B
      rCRS	378	.	C	T	.	clustered_events	        SNP;DP=1612;AF=0.13	SM	chrM.C
      ...
 
    head mutect2.tab
      Run     haplogroup  03%S  03%s  03%I  03%i  05%S  05%s  05%I  05%i  10%S  10%s  10%I  10%i
      chrM.A  A2+(64)     28    35    3     8     28    35    3     8     28    34    3     8
      chrM.B  B2          25    34    2     8     25    34    2     8     25    34    2     8
      chrM.C  C           35    35    4     8     35    35    4     8     35    35    4     8

## LEGEND ##

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
    03%sp			  # homozygous SNPs, 3% heteroplasmy rate, non homopimeric
    ...
    05%
    10%
