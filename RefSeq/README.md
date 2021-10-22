#  SNV Annotation #

## FILES ##

    CDS.bed.gz   : coding regions
    RNR.bed.gz   : rRNA
    TRN.bed.gz   : tRNA
    DLOOP.bed.gz : D-loop
    HV.bed.gz    : hyper variable regions
    HP.bed.gz    : homopolymers

    CADD.vcf.gz  : CADD
    HG.vcf.gz    : haplogroup specific SNV's
    dbSNP.vcf.gz : dbSNP
    NUMT.vcf.gz  : NUMT SNV's
    HS.bed.gz    : hotspots

## SUMMARY ##

    .        #regions  min  q1   q2   q3    max   mean    sum
    CDS      13        207  525  784  1141  1812  876.54  11395
    RNR      2         954  954  1559 1559  1559  1256.50 2513           
    TRN      22        59   68   69   70    75    68.54   1508           
    DLOOP    2         576  576  746  746   746   661.00  1322
    HV       3         136  136  315  359   359   270.00  810
    HP       9         5    9    11   15    22    12.33   111

    .        #SNVs
    CADD     16568*3
    HG       1098
    dbSNP    388
    NUMT     217
    HS       13

## CADD SCORES ##

  .           #      min        q1        q2        q3        max        mean  
  CR(Raw)     49704  -2.198591  0.820158  1.244873  3.607209  11.275703  2.0359	
  CP(PHREDD)  49704  0.001      9.588     11.98     23.2      37         14.9351
