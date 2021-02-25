#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	$opt{i}=1;

	my $result = GetOptions(
		"in=s"	=>	\$opt{in},
		 "i=i"  => 	\$opt{i}
	);
        die "ERROR: $! " if (!$result);

	#######################################################

        my (@samples,%samples);
	open(IN,$opt{in}) or die "ERROR1: $!";
        while(<IN>)
        {
		chomp;
		next if(/^#/ or /^$/);
		my @F=split;
                push @samples,$F[$opt{i}];
		$samples{$F[$opt{i}]}=1;
        }
	close(IN);

	#######################################################
	my ($key,$sample,@keys,%keys,%GT_DP_AF);
	while(<>)
	{
		if(/^##/)
		{
			next if(/^##FILTER/);
			print;
		}
		elsif(/^#/)
		{

			print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
			print "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes\">\n";
		}
		else
		{
			my @F=split;
			$F[2]=".";
			$F[5]=".";
			$F[6]=".";

			my $sample;;
			if($F[7]=~/SM=(.+?);(.+)/)
			{
				$sample=$1;
				$F[7]=$2;
			}
			elsif($F[7]=~/SM=(.+)/)
                        {
                                $sample=$1;
				$F[7]=".";
			}

			defined($samples{$sample}) or die "ERROR3: $_";
			my $key=join "\t",@F[0..7];

			$GT_DP_AF{$key}{$sample}=$F[-1];

			if(!$keys{$key})
			{
				$keys{$key}=1;
				push @keys,$key;
			}
		}
	}

	print join "\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",@samples); print "\n";
	foreach my $key (@keys)
	{
		my @F=();
		my @SF;
		push @F,($key,"GT:DP:AF");

		my $AC=0;
		my $AN=0;
		foreach my $i (1..@samples)
		{
			my $sample=$samples[$i-1];

			if($GT_DP_AF{$key}{$sample})
			{
				push @F,$GT_DP_AF{$key}{$sample} ;
				push @SF,$i;

				my $GT_DP_AF=$GT_DP_AF{$key}{$sample};
				$GT_DP_AF=~/(.+?):/;
				my $GT=$1;
				my @GT=split /[|\/]/,$GT;

				if($GT[0]==1 or $GT[1]==1)
				{
					$AC++;
					$AN+=2;
				}
				else
				{
					#while(@GT and $GT[0] ne "1") { $AN++; shift @GT }
					$AN+=2;
				}
			}
			else { push @F,".:.:." }
		}

		if($F[0]=~/(.+)\t\.$/)       { $F[0]="$1\tAC=$AC;AN=$AN"    }
		elsif($F[0]=~/(.+)\t(\S+)$/) { $F[0]="$1\tAC=$AC;AN=$AN;$2" }
		print join "\t",@F; print "\n";
	}
	exit 0;
}

