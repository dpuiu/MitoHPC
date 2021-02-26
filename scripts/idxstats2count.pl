#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	my %count;

	my $result = GetOptions(
		"sample|Run=s"	=> \$opt{sample}
	);
        die "ERROR: $! " if (!$result);

	while(<>)
	{
		#0	1		2		3
		#chr1	248956422	76763671	382295
		#...
		#chrM	16569		445041		2164

		my @F=split;
		$count{all}+=$F[2]+$F[3];
		$count{mapped}+=$F[2];
		$count{chrM}=$F[2] if($F[0] eq "chrM") 
	}

	print join "\t",("Run","all","mapped","chrM"); print "\n";
	print join "\t",($opt{sample},$count{all},$count{mapped},$count{chrM}); print "\n";

	exit 0;
}
