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
	my (%h0,%h1,$empty,$two);

	# validate input parameters
	my $result = GetOptions(
		"empty=s" =>	\$empty,
		"2"	=> \$two
	);
	die "ERROR: $!" if (!$result);
	###########################################################################

	open(IN,$ARGV[1]) or die "ERROR:$!";
	while(<IN>)
	{
		chomp;
		next if(/^$/);
		my @F=split /\t/;
		my $k1=shift @F;
		$k1.="\t".(shift @F) if(@F and $two);

		$h1{$k1}=join "\t",@F;
	}
	close(IN);

	###########################################################################
	open(IN,$ARGV[0]) or die "ERROR:$!";
	while(<IN>)
	{
		chomp;
		next if(/^$/);
 		my @F=split /\t/;
                my $k0=shift @F;
                $k0.="\t".(shift @F) if(@F and $two);
                $h0{$k0}=join "\t",@F;

		if(defined($h1{$k0}))     { print join "\t",($_,$h1{$k0}); print "\n"; }
		elsif(defined($empty))    { print join "\t",($_,$empty); print "\n"; }
	}
	close(IN);

	if(defined($empty))
	{
		foreach my $k1 (keys %h1)
		{
			if(!defined($h0{$k1})) { print join "\t",($k1,$empty,$h1{$k1}) ; print "\n" }
		}
 	}

	exit 0;
}
