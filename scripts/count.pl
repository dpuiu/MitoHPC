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
	my $j;
	my $min;
	my $max;
	my %count;

        my $result = GetOptions(
		"i=i"	=> \$i,
		"j=i"	=> \$j,
		"min=i"	=> \$min,
		"max=i"	=> \$max
	);
	die "ERROR" if (!$result);


	######################################################
	while(<>)
	{
		next if(/^$/ or /^#/) ;

		chomp;
		#my @F=split /\t/;
		my @F=split;

		next if(@F<=$i);
		next if(defined($j) and @F<=$j);

		if(defined($j)) { $count{$F[$i]}+=$F[$j] }
		else		{ $count{$F[$i]}++;      }
	}

	##########################################################

	foreach my $key (sort {$count{$b}<=>$count{$a}} keys %count)
	{
		next if(defined($min) and $count{$key}<$min);
		next if(defined($max) and $count{$key}>$max);

                print $key,"\t", $count{$key};
                print "\n";
	}

	exit 0;
}

