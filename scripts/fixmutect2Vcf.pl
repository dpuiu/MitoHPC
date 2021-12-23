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
        my $result = GetOptions(
		"file=s"         => \$opt{file},
	);
        die "ERROR: $! " if (!$result);

	############################################################

	while(<>)
	{
		print;
	}
}

