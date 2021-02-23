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

	while(<>)
	{
		chomp;
		next if(/^#/ or /^$/);
		my @F=split;

		print "$F[0]\t$1\n" if($F[0]=~/.+\/(\S+)\./ or $F[0]=~/(\S+)\./)
	}
	exit 0;
}

