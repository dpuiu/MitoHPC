#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters unique lines based on 2 columns

        EXAMPLE:
                 cat I.vcf | uniq2.pl -i 1 -j -1
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %h;
	my ($i,$j)=(0,1);

        my $result = GetOptions(
		"i=i"	=> \$i,
                "j=i"   => \$j
	);
	die "ERROR: $! " if (!$result);

	while(<>)
	{
		chomp;
		my @F=split /\t/;
		if(/^$/ or /^#/) { print; next}
		if($i>=@F or $j>=@F) { next }

		my $key="$F[$i] $F[$j]\n";
		print $_,"\n" unless($h{$key});

		$h{$key}=1;
	}

	exit 0;
}
