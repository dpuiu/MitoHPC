#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions(
		"sm=s"	=>	\$opt{sm}
	);
        die "ERROR: $! " if (!$result);

	#0	1	2	3	4	5	6	7		8			9
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO		FORMAT			SAMPLE
	#chrMT	16318	.	AG	A	.	PASS	AC=1,1,1;AN=4	GT:DP:HF:CILOW:CIUP:SDP	0/1/0/0:1652:0.185,0.004,0.003:0.167,0.002,0.001:0.204,0.008,0.007:151;154,7;0,0;6
	#chrMT	16318	.	AG	CG	.	PASS	AC=1,1,1;AN=4	GT:DP:HF:CILOW:CIUP:SDP	0/0/1/0:1652:0.185,0.004,0.003:0.167,0.002,0.001:0.204,0.008,0.007:151;154,7;0,0;6
	#chrMT	16318	.	AG	TG	.	PASS	AC=1,1,1;AN=4	GT:DP:HF:CILOW:CIUP:SDP	0/0/0/1:1652:0.185,0.004,0.003:0.167,0.002,0.001:0.204,0.008,0.007:151;154,7;0,0;6

	while(<>)
	{
		if(/^#/)
		{
			print;
			next
		}

		my @F=split;

		$F[9]=~/^(.+?):(.+?):(.+?):/ or die "ERROR $_";
		my @F9=split /:/,$F[9];

		my @GT=split /\//,$F9[0];
		my $DP=$F9[1];
		my @AF=split /,/,$F9[2];
		foreach my $i (0..@GT-1)
		{
			if($GT[$i+1])
			{
				$F[9]="$1:$2:$AF[$i]";
				last;
			}
		}
		$F[7]="SM=$opt{sm}";
		$F[8]="GT:DP:AF";
		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

