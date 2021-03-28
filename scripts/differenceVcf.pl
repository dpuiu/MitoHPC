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
        my %h;

        # validate input parameters
        my $result = GetOptions();

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
		#0	1	2	3	4	5	6	7				8	9
		#chrM	73	.	A	G	.	PASS	HV;GT=0/1;DP=360;AF=0.354	SM	M1-53_S82_L001_sorted

                next if(/^$/ or /^#/);

                my @F=split;
		my $SM="";

                if(@F>=10)
		{
			if($F[8] eq "SM")                                 { $SM=$F[9] }
			elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
		}

                $h{"$F[0] $F[1] $F[3] $F[4] $SM"}=1;
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

                my @F=split;
               	my $SM="";

                if(@F>=10)
               	{
                        if($F[8] eq "SM")                                 { $SM=$F[9] }
                        elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
               	}

                print unless $h{"$F[0] $F[1] $F[3] $F[4] $SM"};
        }

	exit 0;
}


