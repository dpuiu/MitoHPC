#!/usr/bin/perl -w

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
        my ($i,$j)=(0,0);
        my %h;

        # validate input parameters
        my $result = GetOptions(
                "i=i"  => \$i,
                "j=i"  => \$j,
        );

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);

                my @F=split;
                die "ERROR:$_" if(@F<$j);
                $h{$F[$j]}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

                my @F=split;
                die "ERROR:$_" if(@F<$i);
                print unless $h{$F[$i]};
        }

	exit 0;
}


