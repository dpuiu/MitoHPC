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
        my $result = GetOptions(
		"sm"	=>	\$opt{sm},
		"min1=s"	=>	\$opt{m1},
		"Max1=s"  =>      \$opt{M1},
               	"min2=s"  =>      \$opt{m2},
               	"Max2=s"  =>      \$opt{M2},
		"alt"	=>	\$opt{alt}
	);

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

                if(@F>=10 and $opt{sm})
		{
			if($F[8] eq "SM")                                 { $SM=$F[9] }
			elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
			else						  { die "ERROR: $_" }
		}

		my $AF=1;
		$AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);

		next if(defined($opt{m2}) and $AF<$opt{m2});
		next if(defined($opt{M2}) and $AF>$opt{M2});

                $h{"$F[0] $F[1] $F[3] $F[4] $SM"}=1;
		$h{"$F[0] $F[1] $F[4] $F[3] $SM"}=1 if($opt{alt});
        }
	close(IN);
        #last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^$/ or /^#/) { print; next}

                my @F=split;
               	my $SM="";

                if(@F>=10 and $opt{sm})
               	{
                        if($F[8] eq "SM")                                 { $SM=$F[9] }
                        elsif($F[7]=~/SM=(\S+?);/ or $F[7]=~/SM=(\S+?)$/) { $SM=$1 }
			else                                              { die "ERROR: $_" }

               	}

		my $AF=1;
                $AF=$1 if(/AF=(0\.\d+)/ or /.+:(1)$/ or /.+:(0\.\d+)/);

               	next if(defined($opt{m1}) and $AF<$opt{m1});
                next if(defined($opt{M1}) and $AF>$opt{M1});

                print unless $h{"$F[0] $F[1] $F[3] $F[4] $SM"};
        }

	exit 0;
}


