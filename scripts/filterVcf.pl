#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that filters a VCF file

        EXAMPLE:
                 cat I.vcf | filterVcf.pl -sample SRR00000 -source mutect2 -p 0.05
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %opt;
	$opt{percent}=0;

	my $result = GetOptions(
		"sample=s"	=> \$opt{sample},
		"source=s"	=> \$opt{source},
                "percent=s"     => \$opt{percent},
        );
        die "ERROR: $! " if (!$result);

	while(<>)
	{
		chomp;s/\|/\//g;
                my @F=split /\t/;

		if(/^#CHROM/)
		{
		 	$F[9]=$opt{sample} if($opt{sample});
			print join "\t",@F[0..9];
			print "\n";
			next;
		}
               	elsif(/^#/)
		{
			print;
			print "\n";
			next;
		}
		if($F[8]=~/^GT:AF:DP$/ and $F[9]=~/^(.+?):(.+?):(.+?)$/ or $F[8]=~/^GT:AF:DP:/ and $F[9]=~/^(.+?):(.+?):(.+?):/ or $F[8]=~/^GT:AD:AF:DP:/ and  $F[9]=~/^(.+?):.+?:(.+?):(.+?):/)
		{
			my ($GT,$AF,$DP)=($1,$2,$3);
			my @AF=split /,/,$AF;
			if(@AF>1)
			{
				$GT=~/(\d+)/;
				$GT=$1;
				$GT=$GT-- if($GT>1);
				$AF=$AF[$GT];
			}

			$F[7]=($F[4] eq "*" or length($F[3]) ne length($F[4]))?"INDEL":"SNP";
			if($AF<$opt{percent})
			{
				next;
			}
			elsif($AF>=$opt{percent} or $AF>=1-$opt{percent})
			{
				$F[7].=";DP=$DP;AF=$AF";
			}
			else
			{
				$F[7].=";DP=$DP";
			}
		}
		elsif($F[8]=~/^GT:DP$/ and $F[9]=~/^(.+?):(.+?)$/ or $F[8]=~/^GT:DP:/ and $F[9]=~/^(.+?):(.+?):/)
		{
			my ($GT,$DP)=($1,$2);

			$F[7]=($F[4] eq "*" or length($F[3]) ne length($F[4]))?"INDEL":"SNP";
			$F[7].=";DP=$DP";
		}
		elsif($F[7]=~/(.+);AF=(\d\.\d+)(.*)/ or $F[7]=~/(.+);AF=(\d)(.*)/)
                {
			my $AF=$2;
			if($AF<$opt{percent})
			{
				next;
                        }
			elsif($AF>1-$opt{percent})
			{
				$F[7]="$1$3";
			}
		}

		@F[8,9]=("SM",$opt{sample}) if($opt{sample});

		my @F3=split //,$F[3];
		my @F4=split //,$F[4];
		if(scalar(@F3)==scalar(@F4) and scalar(@F3)>1)
		{
			foreach my $i (0..@F3-1)
			{
				next if($F3[$i] eq $F4[$i]);
				print join "\t",($F[0],$F[1]+$i,$F[2],$F3[$i],$F4[$i],@F[5..9]); 
				print "\n";
			}
		}
		else
		{
			print join "\t",@F[0..9];print "\n";
		}
	}
	exit 0;
}

