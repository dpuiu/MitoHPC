paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl dbSNP.chrM.vcf - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > dbSNP.RSRS.vcf
paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl CADD.chrM.vcf  - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > CADD.RSRS.vcf
paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl NUMT.chrM.vcf  - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > NUMT.RSRS.vcf

~/sw/bin/minimap2 ../examples1/RefSeq/chrM.30.fa NUMT51.fa -a | samtools view -b | bedtools bamtobed -ed| ../scripts/bed2bed.pl -ed | sort -k5,5nr | uniq.pl -i 3 | ~/bin/join.pl - NUMT51.fa.fai -i 3| cut -f1,2,3,4,5,6,7 | column -t | head
chrM.L3  4630   16469  KM281521.1                 11776  +  12416
chrM.R   0      10933  KM281524.1                 10911  +  16759
chrM.HV  0      9527   KM281533.1                 9525   +  14317
chrM.L3  3913   9755   chr1:628834-635104         5760   +  6271
chrM.HV  0      3937   KM281532.1                 3934   +  7763
chrM.A   12068  13735  KM281526.1                 1642   +  2027
chrM.I   7078   8495   KM281523.1                 1369   +  1870
chrM.L1  869    1533   chr5:80651184-80651847     640    -  664
chrM.A   11610  12262  chr5:134926533-134927184   628    -  652
chrM.L3  16023  16569  KM281514.1                 533    -  626

