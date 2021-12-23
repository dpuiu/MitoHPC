paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl dbSNP.chrM.vcf - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > dbSNP.RSRS.vcf
paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl CADD.chrM.vcf  - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > CADD.RSRS.vcf
paste chrM.1  RSRS.1  | nl |  ~/bin/join.pl NUMT.chrM.vcf  - -i 1 | perl -ane 'next if($F[4] eq $F[-1]); print' | sed 's|chrM|RSRS|'  | cut -f1,2,3,4,5,6,7,8 > NUMT.RSRS.vcf
