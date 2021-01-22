#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that adds the homopolimer(HP) label to certain VCF file SNPs

        EXAMPLE:
                 cat I.vcf | labelVcf.pl 
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
	my @begin=(296,493,511,540,564,954,5888,8269,16178) ; #(300,511,564,2456,4604,5888,6691,8269,8490,9477,11031,12417,13230,14503,16178);
	my @end=(318,502,524,545,574,965,5895,8288,16193);    #(317,524,573,2463,4611,5895,6698,8288,8502,9484,11038,12425,13237,14510,16193);

	my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	while(<>)
	{
		if(/^#/) { print; next }

		chomp;
		s/HP;//;
		my @F=split /\t/;

		foreach my $i (0..@begin-1)
		{
			if($begin[$i]<$F[1] and $F[1]<=$end[$i])
			{
				$F[7].=";HP";
				#last;
			}
		}

		print join "\t",@F; print "\n";
	}

	exit 0;
}

