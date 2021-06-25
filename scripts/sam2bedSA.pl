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
	# define variables
	my %options;
	my %len;

	# validate input parameters
	my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	####################################

	while(<>)
	{
		my @F=split;
                next unless(@F);

		if(/^\@SQ\s+SN:(\S+)\s+LN:(\S+)/)
		{
			$len{$1}=$2;
			next;
		}
		elsif(/^\@/ or @F<11)
                {
                        next;
                }

		next if($F[1]!~/^\d+$/);
		next if($F[1] & 0x4);
		next if($options{mated} and !($F[1] & 0x2) );
		next if($F[2] eq "*");
		next if($F[1] and $F[1] & 0x100) ;  #secondary alignment

		#################################################

		#SA:Z:chrM,6577,-,101S50M,60,0;
		if(/\tSA:Z:($F[2],\S+)/)
		{
			my $sa=$1;

			#################################################
			my ($ref,$begin,$end,$CIGAR,$cigar,$qry,$strand,$score);

	                $qry=$F[0];
	                $qry=$qry."/1" if($F[1] & 0x40);
        	        $qry=$qry."/2" if($F[1] & 0x80);
                	$strand=($F[1] & 0x10 )?"-":"+";
			$ref=$F[2];

	                $CIGAR=$cigar=$F[5];
        	        $begin=$end=$F[3]-1 ;
                	while($cigar and $cigar=~/(\d+)(\w)(.*)/)
                	{
                        	$end+=$1 if($2 eq "M" or $2 eq "D" or $2 eq "N") ;
                        	$cigar=$3;
               	 	}
                	$score=$end-$begin;

			#################################################

			$sa=~/(\w+),(\d+),(.),(\w+),(\d+),(\d+);(.*)/;
			my ($saref,$sabegin,$saend,$saCIGAR,$sacigar,$saqry,$sastrand,$sascore);

			$saref=$1;
			$sastrand=$3;
			$sabegin=$saend=$2-1;
			$saCIGAR=$sacigar=$4;

		        while($sacigar and $sacigar=~/(\d+)(\w)(.*)/)
                	{
                        	$saend+=$1 if($2 eq "M" or $2 eq "D" or $2 eq "N") ;
                       		$sacigar=$3;
                	}
			$sascore=$saend-$sabegin;

			if($begin<$sabegin)
			{
				print join "\t",($ref,$begin,$end,$qry,$score,$strand,$CIGAR);print "\t";
				print join "\t",($saref,$sabegin,$saend,$qry,$sascore,$sastrand,$saCIGAR);print "\n";
			}
			else
                        {
				print join "\t",($saref,$sabegin,$saend,$qry,$sascore,$sastrand,$saCIGAR);print "\t";
                                print join "\t",($ref,$begin,$end,$qry,$score,$strand,$CIGAR);print "\n";
                         }
		}
	}
	exit 0;
}

