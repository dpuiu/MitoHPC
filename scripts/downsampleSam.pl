#!/usr/bin/perl -w

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

	my (%h,%k,%l);

	# read SAM file
	while(<>)
	{
		if(/^#/ or /^$/) {}
		elsif(/^\@/) 
		{ 
			if(/^\@SQ\tSN:(\S+)\s+LN:(\d+)/)
			{
				$l{$1}=$2;
			}
			print ;
		}
		else
		{
			my @F=split;

			die if(@F<8);
			$F[6]=$F[2] if($F[6] eq "=");
			
			my $k=($F[3]<=$F[7])?1:0;

                        if($hg38)
                        {                                        
				if($F[2] eq "chrM"  and $F[3]==1 and $F[5]=~/^(\d+)[SH]/)        {               $F[3]=$l{$F[2]}-$1 }
				if($F[2] eq "chr1"  and $F[3]>=629084       and $F[3]<=634672)   { $F[2]="chrM"; $F[3]-=629084-3914-1 }
				if($F[2] eq "chr17" and $F[3]>=22521208     and $F[3]<=22521400) { $F[2]="chrM"; $F[3]-=22521208-16377-1 }
                                if($F[2] eq "chr17" and $F[3]>=22521401     and $F[3]<=22521639) { $F[2]="chrM"; $F[3]-=22521401-1-1  }			 


				if($F[6] eq "chrM"  and $F[7]==1 and /\tMC:Z:(\d+)[SH]/)	 {               $F[7]=$l{$F[6]}-$1 }
                                if($F[6] eq "chr1"  and $F[7]>=629084       and $F[7]<=634672)   { $F[6]="chrM"; $F[7]-=629084-3914-1 }
				if($F[6] eq "chr17" and $F[7]>=22521208     and $F[7]<=22521400) { $F[6]="chrM"; $F[7]-=22521208-16377-1 }
                                if($F[6] eq "chr17" and $F[7]>=22521401     and $F[7]<=22521639) { $F[6]="chrM"; $F[7]-=22521401-1-1 }			
                        }

			if($k{$F[0]}) 
			{ 
				print ;
			}							   
			elsif($k)
			{
				next if($hg38 and !($F[2] eq "chrM" and $F[6] eq "chrM"));

				$F[3]=$l{$F[2]}-$1 if($F[3]==1 and $F[5]=~/^(\d+)[SH]/);

				#check read count
				next if($h{"$F[2] $F[3]"} and $h{"$F[2] $F[3]"}>=$max);
                                next if($h{"$F[6] $F[7]"} and $h{"$F[6] $F[7]"}>=$max);

				#increment counts
				$h{"$F[2] $F[3]"}++;
				$h{"$F[6] $F[7]"}++;

				#save read name
				$k{$F[0]}=1;

				print;
			}
		}
	}

	exit 0;
}
