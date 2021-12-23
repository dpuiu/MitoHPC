#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	#0	1	2	3	4	5	6	7			8		9
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO			FORMAT	
	#chrM	64	.	C	T	.	PASS	SNP;DP=653;AF=0.038;HV	SM		SRR0000000
	#chrM   15	.	C	T	.       PASS    SM=SRR0000000;SNP;NUMT  GT:DP:AF        0/1:97:1

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
			s/INFO/FORMAT/ if(/^##INFO=<ID=DP,/ or /##INFO=<ID=AF,/);
			s/FORMAT/INFO/ if(/^##FORMAT=<ID=SM,/); 

			print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" if(/^##FORMAT=<ID=DP,/);
			print;
			next;
		}
		
		my @F=split;
		if($F[8] eq "SM")
		{
			my ($GT,$DP,$AF,$ANNOTATION);

			#GT=0/1;DP=60;AF=0.051;DLOOP:CR=1.056772;CP=10.98

			$ANNOTATION="SM=$F[-1]";
			($F[7]=~/(.*)GT=(.+?);DP=(\d+);AF=(\d+\.\d+)(.*)/ or $F[7]=~/(.*)GT=(.+?);DP=(\d+);AF=(\d+)(.*)/) or die "ERROR $_";
			{
				$GT=$2;
				$DP=$3;
				$AF=$4;

				$ANNOTATION.=";$1" if($1); 
				$ANNOTATION.="$5"  if($5);
				$ANNOTATION=~s/;;/;/;
			}

			$F[7]=$ANNOTATION;
			$F[8]="GT:DP:AF";
			$F[9]="$GT:$DP:$AF";
		}

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

