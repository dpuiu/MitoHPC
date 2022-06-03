#!/usr/bin/env perl 

use strict;
use warnings;
use Getopt::Long;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my (%opt,%h);
	my($i,$j)=(0,0);	#-1,0

        # validate input parameters
        my $result = GetOptions(
		"i=i"	=>\$i,
		"j=i"	=>\$j,
	);

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^\@/ or /^#/);

                my @F=split;
                $h{"$F[$j]"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^\@/ or /^#/) { print; next}

                my @F=split;
		print unless($h{$F[$i]});
        }
	close(IN);

	exit 0;
}


