#!/usr/bin/perl -w

use strict;
use Getopt::Long;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my $i;
	my $max=1;

        my $result = GetOptions(
                "i=i"   => \$i,
		"max=i"	=> \$max
        );
	die "ERROR: $? " if (!$result);

	my %h;
	while(<>)
	{
		if(defined($i))
		{
			if(/^#/ or /^$/) { print }
			elsif(/^@/) { print }
			else
			{
				my @F=split;
				die if(!defined($F[$i]));

				$h{$F[$i]}++;
				print if($h{$F[$i]}<=$max);
			}
		}
		else
		{
			$h{$_}++;
                        print if($h{$_}<=$max);
		}
	}

	exit 0;
}

