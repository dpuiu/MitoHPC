#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that counts ...

        EXAMPLE:
                 cat I.vcf | reduceVcf.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my %opt;
        my $result = GetOptions();
        die "ERROR: $! " if (!$result);

	############################################################
	#0	1	2	3	4	5	6	7		8	9
	#chrM	3	.	T	C	.	PASS	SNP;DP=59;HV	SM	chrM.L0
	#=>
	#chr	pos	ref	alt	filters	AC_hom	AC_het	AF_hom		AF_het		AN	max_observed_heteroplasmy
	#chrM	3	T	C	PASS	19	1	3.3667646E-4	1.7719814E-5	56434	9.9700e-01

	my %h;
	while(<>)
	{
		if(/^#/)
		{
			print if(/^##reference=/ or /^##contig/);
			print "$1\n" if(/(^#CHROM.+INFO)/);
			next;
		}

		chomp;
		my @F=split /\t/;
		$F[6]=".";

		my $key="$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]";

		$F[7]=~/(.+);DP=(\d+)/ or die "ERROR: $_";

		$h{$key}{"INFO"}=$1;
		$h{$key}{"DPA"}+=$2;

		if($F[7]=~/AF=(0\.\d+)(.*)/)
		{
			$h{$key}{"AFC"}++;
			$h{$key}{"AFM"}=$1   if(!$h{$key}{"AFM"} or $1>$h{$key}{"AFM"});
			$h{$key}{"INFO2"}=$2 if($2);
		}
		elsif($F[7]=~/DP=\d+(.*)/)
		{
			$h{$key}{"AC"}++;
			$h{$key}{"INFO2"}=$1 if($1);
		}
	}

	foreach my $key (keys %h)
	{
		foreach("AC","AFC")
		{
			$h{$key}{$_}=0 unless($h{$key}{$_});
		}
		$h{$key}{"DPA"}=int($h{$key}{"DPA"}/($h{$key}{"AC"}+$h{$key}{"AFC"})+.5);

		foreach("DPA","AC","AFC","AFM")
		{
			$h{$key}{"INFO"}.=";$_=$h{$key}{$_}" if($h{$key}{$_});
		}
		$h{$key}{"INFO"}.=$h{$key}{"INFO2"} if($h{$key}{"INFO2"});
		print $key,"\t",$h{$key}{"INFO"},"\n";
	}

}

