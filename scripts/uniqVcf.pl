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

		my $SM="";
		my $AF=1;
		if(@F>8)
		{
			if($F[8] eq "SM") { $SM=$F[9]}
			elsif($F[7]=~/SM=(.+?);/ or $F[7]=~/SM=(.+)$/) { $SM=$1}

			$AF=$1 if($F[7]=~/AF=(\S+?);/ or $F[7]=~/AF=(\S+)$/);
		}
	
		my $key=join "\t",(@F[0..4],$SM);

		if(!$max{$key} or $AF>$max{$key})
		{
			$max{$key}=$AF;
			$line{$key}=$_;
		}
	}

	foreach my $key ( keys %max )
	{
		print $line{$key};
	}
}

