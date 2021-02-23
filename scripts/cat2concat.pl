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
	#chrM	64	.	C	T	.	PASS	SNP;DP=653;AF=0.038;HV	SM		SRR6771205
	#chrM   15	.	C	T	.       PASS    SM=SRR8013160;SNP;NUMT  GT:DP:AF        0/1:97:1

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
			s/INFO/FORMAT/ if(/^##INFO=<ID=DP,/ or /##INFO=<ID=AF,/);
			s/FORMAT/INFO/ if(/^##FORMAT=<ID=SM,/); 

			print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" if(/^##FORMAT=<ID=DP,/);
			print;
			next
		}
		my @F=split;
		my ($POS,$SM,$DP,$AF,$GT,$ANNOTATION,$INDEL)=(@F[1,-1],0,1,"0/","",0);

		$GT{$POS}{$SM}++;

		foreach (1..$GT{$POS}{$SM}-1) { $GT.="0/" }
		$GT.="1";

		$INDEL=1 if($F[7]=~/INDEL/);
		if($F[7]=~/DP=(\d+);AF=(\d+\.\d+)(.*)/ or $F[7]=~/DP=(\d+);AF=(\d+)(.*)/)
		{
			$DP=$1;
			$AF=$2;
			$ANNOTATION=$3;
		}
		elsif($F[7]=~/DP=(\d+)(.*)/)
		{
			$DP=$1;
			$ANNOTATION=$2
		}
		else
		{
			die "ERROR:$_"
		}

		$F[7]="SM=$SM";
		$F[7].=";INDEL" if($INDEL);
		$F[7].=$ANNOTATION if($ANNOTATION);

		$F[8]="GT:DP:AF";
		$F[9]="$GT:$DP:$AF";

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

