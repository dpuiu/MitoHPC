#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that recalculates snp positions

        EXAMPLE:
		cat I.mutserve.vcf | fixmutserveVcf.pl -file rCRS.fa
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
		"ref=s"		=> \$opt{ref},
 	        "file=s" 	=> \$opt{file},
 	);
	die "ERROR: $! " if (!$result);
	die "ERROR: $! " if (!defined($opt{file}));

	my %diff;
	open(IN,$opt{file}) or die "ERROR:$!";
        while(<IN>)
        {
		if(/INDEL;/)
		{
			my @F=split;
			$diff{$F[1]}=length($F[4])-length($F[3]);
		}
        }
	close(IN);

	if($opt{ref} eq "rCRS") { $diff{3107}=-1 }

	my $diff=0;
	my @pos=(sort {$a<=>$b} keys %diff);
	foreach my $pos (@pos)
	{
		$diff+=$diff{$pos};
		$diff{$pos}=$diff;
	}

	###################################################
	##reference=file://out/chrM.L0/chrM.L0.mutect2.fa>
	##contig=<ID=chrM.L0,length=16567>

	while(<>)
	{
		if(/^$/)
		{
			next;
		}
		elsif(/^##reference/ or /^##contig/)
		{
		}
		elsif(/^#/)
		{
			print;
			next;
		}
		else
		{
			chomp;
			my @F=split /\t/;
			my $diff=0;
			my $i=0;
			while($i<@pos and $F[1]>$pos[$i])
			{
				$diff=$diff{$pos[$i]};
				$i++
			}

			$F[0]=$opt{ref};
			$F[1]-=$diff;
			print join "\t",@F;
			print "\n";
		}
	}
}

