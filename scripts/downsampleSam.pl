#!/usr/bin/perl -w

use strict;
use Getopt::Long;

###############################################################################
#
# Main program
#
# Programs which reduces 
###############################################################################

MAIN:
{
	my ($max,$and,$or)=(1,0,0);
	

        my $result = GetOptions(
		"max=i"	=> \$max,	# max start position count (Ex:10)
		"and" 	=> \$and,	# if both   mates alignment start position count is > $max, the mates are  filtered out
		"or" 	=> \$or,	# if either mates alignment start position cvg is > $max, the mates are filtered out
        );
	die "ERROR: $? " if (!$result);

	my (%h,%k);
	while(<>)
	{
		if(/^#/ or /^$/ or /^\@/) { print }
		else
		{
			my @F=split;

			die if(@F<8);
			$F[3]-=$1 if($F[5]=~/^(\d+)S/); 

			if($k{$F[0]}) { print }
			elsif($F[6] eq "=" and $F[3] and $F[7] and $F[3]<=$F[7])
			{
				next if($h{"$F[2] $F[3] $F[7]"} and $h{"$F[2] $F[3] $F[7]"}>=1) ; 
				$h{"$F[2] $F[3] $F[7]"}++;

				next if($and and ($h{"$F[2] $F[3]"} and $h{"$F[2] $F[3]"}>=$max and $h{"$F[2] $F[7]"} and $h{"$F[2] $F[7]"}>=$max));
				next if($or  and ($h{"$F[2] $F[3]"} and $h{"$F[2] $F[3]"}>=$max or $h{"$F[2] $F[7]"} and $h{"$F[2] $F[7]"}>=$max ));
				$h{"$F[2] $F[3]"}++;
				$h{"$F[2] $F[7]"}++;
				$k{$F[0]}=1;
				print;
			}
		}
	}

	exit 0;
}
