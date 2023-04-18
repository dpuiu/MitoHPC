#!/usr/bin/env perl

use strict;
use Getopt::Long;

###############################################################################
#
# Program which downsamples a SAM file based on read alignment start positions
#   max (default 1) alignments are allowed to start at a certain position;
#    exceptions allowed for the mates of the reads already included
#
# Example:
#
#   samtools view -h run.bam chrM | perl downsampleSam.pl -m 10 | satools view -b > run.downsample.bam
#   samtools index run.downsample.bam

###############################################################################

MAIN:
{
	my $max=10;
	my $hg38=1;

	my $result = GetOptions(
		"max=i"	=> 	\$max,	
		"hg38"	=>	\$hg38
        );
	die "ERROR: $? " if (!$result);

	my (%count,%keep,%skip,%len);

	# read SAM file
	while(<>)
	{
		if(/^#/ or /^$/) {}
		elsif(/^\@/) 
		{ 
			if(/^\@SQ\tSN:(\S+)\s+LN:(\d+)/)
			{
				$len{$1}=$2;
			}
			print ;
		}
		else
		{
			my @F=split;

			die if(@F<8);
			$F[6]=$F[2] if($F[6] eq "=");
			
			my ($F3,$F7)=@F[3,7];
		
                        if($hg38)
                        {                                        
				if($F[2] eq "chr1"  and $F[3]>=629084       and $F[3]<=634672)   { $F[2]="chrM"; $F[3]-=629084-3914-1 }
				if($F[2] eq "chr17" and $F[3]>=22521208     and $F[3]<=22521400) { $F[2]="chrM"; $F[3]-=22521208-16377-1 }
                                if($F[2] eq "chr17" and $F[3]>=22521401     and $F[3]<=22521639) { $F[2]="chrM"; $F[3]-=22521401-1-1  }			 
				if($F[5]=~/^(\d+)[SH]/)                                          {               $F[3]-=$1; $F[3]+=$len{$F[2]} if($F[3]<=0); }

                                if($F[6] eq "chr1"  and $F[7]>=629084       and $F[7]<=634672)   { $F[6]="chrM"; $F[7]-=629084-3914-1 }
				if($F[6] eq "chr17" and $F[7]>=22521208     and $F[7]<=22521400) { $F[6]="chrM"; $F[7]-=22521208-16377-1 }
                                if($F[6] eq "chr17" and $F[7]>=22521401     and $F[7]<=22521639) { $F[6]="chrM"; $F[7]-=22521401-1-1 }			
				if(/\tMC:Z:(\d+)[SH]/)                                           {               $F[7]-=$1; $F[7]+=$len{$F[2]} if($F[7]<=0); }
                        }

			next if($hg38 and !($F[2] eq "chrM" and $F[6] eq "chrM"));

			if($keep{$F[0]}) 
			{ 
				print ; 
			}							   
			elsif(!$skip{$F[0]})
			{
				#check read count
				#if($count{"$F[2] $F[3]"} and $count{"$F[2] $F[3]"}>=$max and $count{"$F[6] $F[7]"} and $count{"$F[6] $F[7]"}>=$max) 
				if($count{"$F[2] $F[3]"} and $count{"$F[2] $F[3]"}>=$max or $count{"$F[6] $F[7]"} and $count{"$F[6] $F[7]"}>=$max)
				{ 
					$skip{$F[0]} = 1;
				}
				else
				{
					#increment counts
					$count{"$F[2] $F[3]"}++;
					$count{"$F[6] $F[7]"}++;

					#save read name
					$keep{$F[0]}=1;

					print;
				}
			}
		}
	}

	exit 0;
}
