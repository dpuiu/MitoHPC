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

        # validate input parameters
        my $result = GetOptions();

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);

                my @F=split;
                $h{"$F[0]"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^@/) { print; next}

                my @F=split;
		print if($h{$F[0]});
        }
	close(IN);

	exit 0;
}


