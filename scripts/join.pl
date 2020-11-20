#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that joins 2 tables by the first column

        EXAMPLE:
                 join I.tab J.tab -empty 0
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %h;
	my $empty;

	# validate input parameters
	my $result = GetOptions(
		"empty=s" =>	\$empty,
	);
	die "ERROR: $!" if (!$result);
	###########################################################################

	open(IN,$ARGV[1]) or die "ERROR:$!";
	while(<IN>)
	{
		next if(/^$/ or /^#/);
		chomp;
		my @F=split /\t/;
		my $F=shift @F;
		$h{$F}=join "\t",@F;
	}
	close(IN);

	###########################################################################
	open(IN,$ARGV[0]) or die "ERROR:$!";
	while(<IN>)
	{
		next if(/^$/ or /^#/);
		chomp;
		my @F=split /\t/;
		my $F=$F[0];

		if(defined($h{$F}))     { print join "\t",(@F,$h{$F}); print "\n"; }
		elsif(defined($empty))  { print join "\t",(@F,$empty); print "\n"; }
	}
	close(IN);
	exit 0;
}
