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
        my (%opt,%id,%info,%info_tag);

        # validate input parameters
        my $result = GetOptions();

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,$ARGV[1]) or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {
		#0	1	2	3	4	5	6	7				
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		#chrM	55	.	T	C	.	.	HG=H

                if(/^##INFO=<ID=(\w+)/) {$info_tag{$1}=$_}
		if(/^#/) { next }

                my @F=split;
		my $key="$F[0] $F[1] $F[3] $F[4]";

		$id{$key}=$F[2]   if($F[2] and $F[2] ne ".");
                $info{$key}=$F[7] if($F[7] and $F[7] ne ".");
        }
	close(IN);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
		if(/^##INFO=<ID=(\w+)/) 
		{ 
			delete($info_tag{$1})
		}
		elsif(/^#CHROM/)
		{
			foreach (keys %info_tag) { print $info_tag{$_} }
		}
                if(/^#/) { print ; next}

                my @F=split;
		my $key="$F[0] $F[1] $F[3] $F[4]";

		if($id{$key} and $F[2] and $F[2]!~$id{$key})
               	{
                       	$F[2]=($F[2] eq ".")?"$id{$key}":"$F[2]:$id{$key}";
               	}
		elsif($info{$key} and $F[7] and $F[7]!~$info{$key})
		{
			$F[7]=($F[7] eq ".")?"$info{$key}":"$F[7]:$info{$key}";
		}
		print join "\t",@F; print "\n";
        }

	exit 0;
}


