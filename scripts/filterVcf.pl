#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters a VCF file

        EXAMPLE:
                 cat I.vcf | filterVcf.pl -sample SRR00000 -source mutect2 -p 0.05
~;

MAIN:
{
	# define variables
	my %opt;
	$opt{percent}="0.0";

	my $result = GetOptions(
                "percent=s" 	=> \$opt{percent},
		"sample=s"	=> \$opt{sample},
		"source=s"	=> \$opt{source},
		"header=s"	=> \$opt{header}
        );
        die "ERROR: $! " if (!$result);
	$opt{percent}=~/^0.[01234]\d*$/ or die "ERROR:percent must be >=0 and <0.5";

	if($opt{header})
	{
		open(IN,$opt{header}) or die "ERROR: $!";
		while(<IN>)
		{
			print;
		}
		close(IN)
	}

	while(<>)
	{
		if(/^#/)
		{
			if($opt{header}) { next }
			else		 { print; next}
		}

		chomp;
                my @F=split /\t/;
		die "NORMALIZATION ERROR:$_\n" if($F[4]=~/,/);
		#$F[6]="." unless($F[6] eq "PASS");

		if($F[7]!~/SNP/ and $F[7]!~/INDEL/)
		{
			$F[7]=($F[4] eq "*" or length($F[3]) ne length($F[4]))?"INDEL":"SNP";
		}

		if($F[7]!~/SM=/ and defined($opt{sample}))
		{
			$F[7]="SM=$opt{sample};$F[7]";
		}

		my @F8=split /:/,$F[8];
		my @F9=split /:/,$F[9];
		my %h;
		foreach my $i (0..@F8-1)
		{
			$h{$F8[$i]}=$F9[$i];
		}
		$h{AF}=1 unless($h{AF});
		$h{AF}=$1 if($h{AF}=~/(\d\.\d+)/);

		if(!$h{AF} or $h{AF}>1-$opt{percent}) 	{ $h{AF}=1; }
		elsif($h{AF}<$opt{percent})		{ next;     }

		$F[8]="GT:DP:AF";
		$F[9]="$h{GT}:$h{DP}:$h{AF}";

		if(length($F[3])>1 and length($F[3])==length($F[4]))
		{
			my @F3=split //,$F[3];
			my @F4=split //,$F[4];

			foreach my $i (0..@F3-1)
			{
				if($F3[$i] ne $F4[$i])
				{
					print join "\t",($F[0],$F[1]+$i,$F[2],$F3[$i],$F4[$i],@F[5..9]);
					print "\n";
				}
			}
		}
		else
		{
			print join "\t",@F[0..9];
			print "\n";
		}

	}
	exit 0;
}

