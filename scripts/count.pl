#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that counts the values in a specific column (default 0)

        EXAMPLE:
                 cat I.vcf | count.pl -i 1
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my $i=0;
	my %count;

        my $result = GetOptions(
		"i=i"		=> 	\$i,
	);
	die "ERROR" if (!$result);


	######################################################
	while(<>)
	{
		next if(/^$/ or /^#/) ;

		chomp;
		my @F=split /\t/;

		next if(@F<=$i);
		$count{$F[$i]}++;
	}

	##########################################################

	foreach my $key (sort {$count{$b}<=>$count{$a}} keys %count)
	{
                print $key,"\t", $count{$key};
                print "\n";
	}

	exit 0;
}

