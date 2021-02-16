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
			#$F[7]=".";
			if($F[7]=~/(.+);SM=.+?;(.+)/)  { $F[7]="$1;$2" } 
			elsif($F[7]=~/^SM=.+?;(.+)/)   { $F[7]="$1"    } 
			elsif($F[7]=~/(.+);SM=/)       { $F[7]="$1"    }
			else			       { $F[7]="."     }
		}

		print join "\t",@F; 
		print "\n";
	}

}

