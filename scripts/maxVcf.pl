#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that computes each SNP major allele 

        EXAMPLE:
                 cat I.vcf | maxVcf.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my (%max,%line,%sum);
	while(<>)
	{
		if(/^#/) { print; next; }

		chomp;
		my @F=split /\t/;

		my $POS=$F[1];

		my $AF=1;
		$AF=$1 if($F[7]=~/AF=(\S+?);/ or $F[7]=~/AF=(\S+)$/);

		if(!$max{$POS} or $AF>$max{$POS})
		{
			$max{$POS}=$AF;
			$line{$POS}=$_;
		}

		$sum{$POS}+=$AF;
	}

	foreach my $POS ( sort {$a <=> $b} keys %max )
	{
		if($max{$POS}>1-$sum{$POS})
		{
			print $line{$POS},"\n";
		}
	}
}

