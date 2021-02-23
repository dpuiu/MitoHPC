#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	my $result = GetOptions(
	);
        die "ERROR0: $! " if (!$result);

	#######################################################

	while(<>)
	{
		if(/^##/)
		{
			next if(/^##FILTER/ or /^##FORMAT/);
			print;
		}
		else
		{
			my @F=split;
			print join "\t",@F[0..7]; print "\n";
		}
	}
	exit 0;
}

