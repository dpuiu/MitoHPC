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
	$opt{m1}=0; $opt{m2}=0; $opt{M1}=1; $opt{M2}=1;

        # validate input parameters
        my $result = GetOptions(
                "min1=s"  =>	  \$opt{m1},
                "Max1=s"  =>	  \$opt{M1},
                "min2=s"  =>	  \$opt{m2},
                "Max2=s"  =>	  \$opt{M2},
	);

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
		#0	1	2	3	4	5	6	7				8	9
		#chrM	73	.	A	G	.	PASS	HV;GT=0/1;DP=360;AF=0.354	SM	M1-53_S82_L001_sorted

                next if(/^$/ or /^#/);

		my $AF=1;
		$AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);
                next if($AF<$opt{m2});
                next if($AF>$opt{M2});

                my @F=split;
                push(@{$h{"$F[0] $F[1] $F[3] $F[4]"}},$_);
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

               	my $AF=1;
                $AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);
                next if($AF<$opt{m1});
                next if($AF>$opt{M1});

                my @F=split;
                foreach my $h (@{$h{"$F[0] $F[1] $F[3] $F[4]"}})
		{
			
			next if($_ eq $h);
			print join "\t",(@F,$h);
		}
        }

	exit 0;
}


