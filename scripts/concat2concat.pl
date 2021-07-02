#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	#0	1	2	3	4	5	6	7			8		9
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO			FORMAT		SAMPLE
	#chrM   15	.	C	T	.       PASS    SM=SRR0000000;SNP;NUMT  GT:AF:DP        0/1:0.99:97
	#=>
	#chrM   15      .       C       T       .       PASS    SM=SRR0000000;SNP;NUMT  GT:DP:AF        0/1:97:0.99

	#0/1:0.798,0.202:104
	#chrM	11251	.	A	T	.	PASS	SM=chrM.T	GT:AF:DP	0/0:0.798,0.202:104
	#chrM	11251	.	A	G	.	PASS	SM=chrM.T	GT:AF:DP	1/0:0.798,0.202:104


	#											0/1/2/3:1790:0.004,0.991,0.005:

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
			print;
			next
		}

		my @F=split;
		my @F8=split /:/,$F[8];
		my @F9=split /:/,$F[9];
		my (%F8);

		foreach my $i (0..@F8-1)
		{
			$F8{$F8[$i]}=$i;
		}

		die "ERROR $_" unless(defined($F8{GT}) and defined($F8{DP}) and defined($F8{AF}));
		$F[8]="GT:DP:AF";

		if($F9[$F8{AF}]=~/(.+),(.+)/)
		{
			if($F9[$F8{GT}] eq "1/0")    { $F9[$F8{AF}]=$1 }
			elsif($F9[$F8{GT}] eq "0/1") { $F9[$F8{AF}]=$2 }
			else                         { next            }
		}

		$F9[$F8{AF}]=1 if($F9[$F8{AF}]==1);

		$F[9]="$F9[$F8{GT}]:$F9[$F8{DP}]:$F9[$F8{AF}]";

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

