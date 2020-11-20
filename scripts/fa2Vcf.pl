#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that generates 3 VCF headers

        EXAMPLE:
                 fa2Vcf.pl F.fa
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	print "##reference=file://$ARGV[0]>\n";

	open(IN,"$ARGV[0].fai") or die "ERROR:$!";
	while(<IN>)
	{
		chomp;
                my @F=split /\t/;
		print "##contig=<ID=$F[0],length=$F[1]>\n";
	}
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

	exit 0;
}

