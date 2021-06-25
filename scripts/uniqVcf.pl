#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	my (%max,%line);
	while(<>)
	{
		if(/^#/)
		{
			print; next;
		}

		my @F=split;

		#0	1	2	3	4	5	6		7				8	9	
		#chrM	567	.	A	ACCCC	.	multiallelic	INDEL;DP=190;AF=0.237;HP	SM	SRR0000000
		#chrM	567	.	A	ACCC	.	multiallelic	INDEL;DP=190;AF=0.095;HP	SM	SRR0000000

		my $POS="$F[1] $F[9]";

		my $AF=1;
		$AF=$1 if($F[7]=~/AF=(\S+?);/ or $F[7]=~/AF=(\S+)$/);

		if(!$max{$POS} or $AF>$max{$POS})
		{
			$max{$POS}=$AF;
			$line{$POS}=$_;
		}
	}

	foreach my $POS ( keys %max )
	{
		print $line{$POS};
	}
}

