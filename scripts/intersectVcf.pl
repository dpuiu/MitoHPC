#!/usr/bin/env perl -w

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
        my ($i0,$i1)=(1,-1);
        my ($j0,$j1)=(1,-1);
        my %h;

        # validate input parameters
        my $result = GetOptions(
                "i0=i"  => \$i0,
                "i1=i"  => \$i1,

                "j0=i"  => \$j0,
                "j1=i"  => \$j1,
        );

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                next if(/^$/ or /^#/);

                my @F=split;
                die "ERROR:$_" if(@F<2);
                $h{"$F[$j0] $F[$j1]"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

                my @F=split;
                die "ERROR:$_" if(@F<2);
                print if $h{"$F[$i0] $F[$i1]"};
        }

	exit 0;
}


