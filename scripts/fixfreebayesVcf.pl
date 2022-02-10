#!/usr/bin/env perl 
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects mutect2 VCF output

        EXAMPLE:
                 cat I.mutect2.vcf | fixmutect2Vcf.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my %opt;
	my (%h,%AF);
        my $result = GetOptions(
		"file=s"         => \$opt{file},
	);
        die "ERROR: $! " if (!$result);

	############################################################

	while(<>)
	{
		#IN
		#0      1       2       3       4       5       6       7       	8                       9
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO		FORMAT			chrM.A
		#chrM	64	.	C	T	15127.1	.	...;AF=; 	GT:DP:AD:RO:QR:AO:QA:GL	1/1:865:0,853:0:0:853:17060:-1523.72,-256.779,0

		#OUT
                #0      1    2   3    4    5     6                                                         7     8         9
		#CHROM  POS  ID  REF  ALT  QUAL  FILTER                                                    INFO  FORMAT    SAMPLE
		#chrM   28   .   A    T    .     clustered_events;strand_bias;strict_strand;weak_evidence  .     GT:DP:AF  0/1:182:0.026

		my @F=split;

		if(/^#CHROM/)
		{
			$F[9]="SAMPLE";
			print join "\t",@F;print "\n";
			next;	
		}
		elsif(/^#/)
		{
			print;
			next
		}
		elsif($F[3]=~/N/)
		{
			next
		}


                $F[5]=int($F[5]+.5);
                $F[7]=".";

                $F[8]=~/^GT:DP:AD:/ or die "ERROR: $_";
                $F[8]="GT:DP:AF";

                $F[9]=~/^(.+?):(.+?):\d+,(.+?):/ or die "ERROR: $_";
                my $AF=int(1000*$3/$2+.5)/1000;
                $F[9]="$1:$2:$AF";


		my $lenDiff=length($F[3])-length($F[4]);
		if($lenDiff<0)
		{
			$F[3]=substr($F[3],0,1);
			$F[4]=substr($F[4],0,-$lenDiff+1);
			print join "\t",@F; print "\n";
		}
		elsif($lenDiff>0)
                {
                        $F[4]=substr($F[4],0,1);
                        $F[3]=substr($F[3],0,$lenDiff+1);
			print join "\t",@F; print "\n";
                }
		else
		{
			foreach my $i (0..length($F[3])-1)
			{
				if(substr($F[3],$i,1) ne substr($F[4],$i,1))
				{
					print join "\t",($F[0],$F[1]+$i,".",substr($F[3],$i,1),substr($F[4],$i,1),@F[5,6,7,8,9]); print "\n";
				}
			}
		}
	}
}

