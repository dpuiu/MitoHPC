#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters unique lines based on 2 columns

        EXAMPLE:
                 cat I.vcf | uniq.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %h;

        my $result = GetOptions(
	);
	die "ERROR: $! " if (!$result);

	while(<>)
	{
		print $_ unless($h{$_});
		$h{$_}=1;
	}

	exit 0;
}
