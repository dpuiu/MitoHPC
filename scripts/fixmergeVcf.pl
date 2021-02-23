#!/usr/bin/env perl 
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects "bcftools merge" output

        EXAMPLE:
                 cat I.mutect2.vcf | fixmergeVcf.pl
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
		"in=s"         => \$opt{in},
	);
        die "ERROR: $! " if (!$result);
	die "ERROR: $! " if (!defined($opt{in}));
	############################################################

	my @I;
	open(IN,$opt{in}) or die "ERROR: $!";
	while(<IN>)
	{
		#chomp;
		#/^\S+$/ or die "ERROR:$_";
		#push @I,$_;

		/.+\/(\S+)\.mut/ or die "ERROR:$_";
		push @I,$1;
	}
	close(IN);
	############################################################

	while(<>)
	{
		#0	1	2	3	4	5	6	7	8	9	10	...
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SMPLE2

		my @F=split;
		if(/^##/)
		{
			print;
			next;
		}

		if(/^#/)
		{
			@F=(@F[0..8],@I);
		}
		else
		{
			$F[6]=".";

			my ($AC,$AN);
			$AC=$1 if($F[7]=~/AC=(\d+)/);
			$AN=$1 if($F[7]=~/AN=(\d+)/);
			if($AC and $AN) { $F[7]="AC=$AC;AN=$AN" }
			elsif($AC)      { $F[7]="AC=$AC"        }
			elsif($AN)      { $F[7]="AN=$AN"        }
		}

		print join "\t",@F; 
		print "\n";
	}

}

