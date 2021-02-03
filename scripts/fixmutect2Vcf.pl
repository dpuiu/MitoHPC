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

	my @P;
	my $PAF;
	while(<>)
	{
		if(/^#/)
		{
			print;
			next;
		}

		chomp;
		my @F=split /\t/;
		if(@P and $P[1]==$F[1] and $P[4] eq $F[4])
		{
			my $PAF=0;
			$P[9]=~/.+:(.+)/;
			$PAF+=$1;

			$F[9]=~/(.+):(.+)/;
			$PAF+=$2;
			$PAF=1 if($PAF>1);

			$F[9]="$1:$PAF";
		}
		elsif(@P)
		{
			print join "\t",@P; 
			print "\n"
		}

		@P=@F
	}

	if(@P)
	{
		print join "\t",@P; 
		print "\n"
	}
}

