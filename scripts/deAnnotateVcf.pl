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
        my @tags=("INDEL","GT","DP","AF","SM");

        # validate input parameters
        my $result = GetOptions(
	);

	die "ERROR: $! " if (!$result);

        #########################################

        while(<>)
        {
                if(/^##INFO=<ID=(\w+),/)
		{
			my $keep=0;
			next unless $1 ~~ @tags;
			print;
		}
		elsif(/^#/)   { print }
		else
		{
	                my @F=split /\t/;
			$F[7]=";$F[7];";
			my $F7="";
			foreach  (@tags)
			{
				$F7.="$1;" if($F[7]=~/;($_);/ or $F[7]=~/;($_=.+?);/);
			}
			$F7="." unless($F7);
			$F7=~s/;$//;
			$F[7]=$F7;

			print join "\t",@F;
		}
        }

	exit 0;
}


