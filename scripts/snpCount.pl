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
                "in=s" 	=> \$opt{in},
		"i=i"	=> \$opt{i}
        );
        die "ERROR: $! " if (!$result);

	#############################################

        my (@samples,%samples,%snp,%count);
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

	############################################

	while(<>)
	{
		next if(/^#/);
		my @F=split;
		$F[7]=~/SM=(.+?);/ or $F[7]=~/SM=(.+)$/ or die "ERROR:$_";
		my $sample=$1;
		defined($samples{$sample}) or die "ERROR: $_";

		if(!$snp{$sample}{$F[1]})
		{
			if($F[-1]=~/.+:1$/)
			{
				$count{$sample}{H}++ ;
				if($F[7]!~/;INDEL/)  { $count{$sample}{S}++ } else { $count{$sample}{I}++ }
			}
			else
			{
				$count{$sample}{h}++ ;
				if($F[7]!~/;INDEL/)  { $count{$sample}{s}++ } else { $count{$sample}{i}++ }
			}

			if($F[7]=~/;HP/)
			{
	        	        if($F[-1]=~/.+:1$/)
        	        	{
                	       	 	$count{$sample}{Hp}++ ;
                        		if($F[7]!~/;INDEL/)  { $count{$sample}{Sp}++ } else { $count{$sample}{Ip}++ }
	                	}
        	        	else
                		{
	                        	$count{$sample}{hp}++ ;
        	                	if($F[7]!~/;INDEL/)  { $count{$sample}{sp}++ } else { $count{$sample}{ip}++ }
	                	}
			}

			$snp{$sample}{$F[1]}=1;
		}
	}

	#########################################################

	print join "\t",("Run","H","h","S","s","I","i","Hp","hp","Sp","sp","Ip","ip"); print "\n";

	foreach my $i (1..@samples)
        {
		my $sample=$samples[$i-1];
		my @counts=();

		foreach ("H","h","S","s","I","i","Hp","hp","Sp","sp","Ip","ip")
		{
			 push @counts,($count{$sample}{$_})?$count{$sample}{$_}:0;
		}

		print join "\t",($sample,@counts); print "\n";
	}

	################################################################
}

