#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	$opt{i}=0;

	my $result = GetOptions(
		"in=s"	=>	\$opt{in},
		"i=i"  => 	\$opt{i}
	);
        die "ERROR: $! " if (!$result);

	my @C;

	#######################################################
	my ($key,$sample,@keys,%keys,%GT_DP_AF);
	while(<>)
	{
		# #cat variants.vcf | grep -v "^##"| cut -f1,2,3,4,5,6,7,8,9,10,11 | head
		# 0             1       2       3       4       5       6       7       8               9               10
		# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT		chrM.V.bam	chrM.B.bam
		# rCRS		64	.	C	T	.	PASS	.	GT:AF:DP	0		0
		# rCRS		72	.	T	C	.	PASS	.	GT:AF:DP	1:1.00:70	0
		# rCRS		73	.	A	G	.	PASS	.	GT:AF:DP	0		1:0.986:69


		if(/^##/)
		{
			print;
		}
		elsif(/^#CHROM/)
		{
			@C=split;
			print join "\t",(@C[0..8],"SAMPLE");print "\n";
		}
		else
		{
			my @F=split;
			foreach my $i (9..@F-1) 
			{ 
				if($F[$i] and $F[$i]=~/(.+?):/ and $1=~/[1-9]/)
				{
					$F[7]="SM=$C[$i]";
					$F[7].=";INDEL" if(length($F[3])!=length($F[4]) or $F[3] eq "*" or $F[4] eq "*");
					print join "\t",(@F[0..8],$F[$i]);
					print "\n";
				}
			}
		}
	}

	exit 0;
}

