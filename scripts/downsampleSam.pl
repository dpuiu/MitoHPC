#!/usr/bin/perl -w

use strict;
use Getopt::Long;

###############################################################################
#
# Program which downsamples a SAM file based on read alignment start positions
#   max (default 1) alignments are allowed to start at a certain position;
#    exceptions allowed for the mates of the reads already included
###############################################################################

MAIN:
{
	my $max=1;

	my $result = GetOptions(
		"max=i"	=> \$max,	# max start position count (Ex:10)
        );
	die "ERROR: $? " if (!$result);

	my (%h,%k);

	# read SAM file
	while(<>)
	{
		if(/^#/ or /^$/ or /^\@/) { print }
		else
		{
			my @F=split;

			die if(@F<8);
			$F[3]-=$1 if($F[5]=~/^(\d+)S/);

			if($k{$F[0]}) { print }						# keep mate
			elsif($F[6] eq "=" and $F[3] and $F[7] and $F[3]<=$F[7])
			{
				#remove duplicates
				next if($h{"$F[2] $F[3] $F[7]"} and $h{"$F[2] $F[3] $F[7]"}>=1) ;
				$h{"$F[2] $F[3] $F[7]"}++;

				#check read count
				next if($h{"$F[2] $F[3]"} and $h{"$F[2] $F[3]"}>=$max);

				#increment counts
				$h{"$F[2] $F[3]"}++;
				$h{"$F[2] $F[7]"}++;

				#save read name
				$k{$F[0]}=1;
				print;
			}
		}
	}

	exit 0;
}
