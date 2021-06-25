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
	#chrM   15	.	C	T	.       PASS    SM=SRR0000000;SNP;NUMT  GT:DP:AF        0/1:97:1
	#chrM   64      .       C       T       .       PASS    SNP;DP=653;AF=0.038;HV  SM              SRR0000000

	my %GT;
	while(<>)
	{
		if(/^#/)
		{
                        s/FORMAT/INFO/ if(/^##FORMAT=<ID=DP,/ or /##FORMAT=<ID=AF,/);
                        s/INFO/FORMAT/ if(/^##INFO=<ID=SM,/); 

			print "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" if(/^##INFO=<ID=DP,/);
			print;
			next
		}

		my @F=split;
		if($F[8] ne "SM")
		{
		
			my ($SM,$ANNOTATION);
			if(/SM=(\S+?);(\S+)/) { ($SM,$ANNOTATION)=($1,"$2;") }
			elsif(/SM=(\S+)\s+/)  { ($SM,$ANNOTATION)=($1,"") }
			else                  { die "ERROR $_";           }

			my ($GT,$DP,$AF)=split /:/,$F[-1];

			$F[7]=$ANNOTATION."GT=$GT;DP=$DP;AF=$AF";
			$F[8]="SM";
			$F[9]="$SM";
		}

		print join "\t",@F;
		print "\n";
	}
	exit 0;
}

