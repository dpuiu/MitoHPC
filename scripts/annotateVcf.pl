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
        my (%info,%h);

        # validate input parameters
        my $result = GetOptions(
	);

	die "ERROR: $! " if (!$result);

        #########################################

        open(IN,"zcat -f $ARGV[1] |") or die("ERROR: Cannot open input file".$!) ;
        while(<IN>)
        {

		##fileformat=VCFv4.2
		##reference=file:rCRS.fa
		##contig=<ID=rCRS,length=16569>
		##INFO=<ID=NUMT,Number=0,Type=Flag,Description="NUMT">
		#0	1	2	3	4	5	6	7
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		#RSRS	15	.	C	T	.	.	NUMT

                if(/^##INFO/) { $info{$_}=""}
		elsif(/^#/)   {}
		else
		{
	                my @F=split;
			my $key="$F[0] $F[1] $F[3] $F[4]";
                	$h{$key}=$F[7];
		}
        }
	close(IN);
        last unless(%h);

        #########################################

        open(IN,$ARGV[0]) or die("ERROR: Cannot open input file".$!) ;

        while(<IN>)
        {
                if(/^#INFO/)
		{
			print $_;
			delete($info{$_});
		}
		elsif(/^#CHROM/)
		{
			print %info;
			print;
		}
		elsif(/^#/)
		{
			print;
		}
		else
		{
	                my @F=split;

			my $key="$F[0] $F[1] $F[3] $F[4]";
			if($h{$key} and $F[7]!~/$h{$key}/)
			{
				$F[7].=";$h{$key}";
			}
			print join "\t",@F;
			print "\n";
		}
        }

	exit 0;
}


