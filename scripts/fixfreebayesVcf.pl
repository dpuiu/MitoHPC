#!/usr/bin/env perl 
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects mutect2 VCF output

        EXAMPLE:
                 cat I.mutect2.vcf | fixmutect2Vcf.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
        my %opt;
	my (%h,%AF);
        my $result = GetOptions(
		"file=s"         => \$opt{file},
	);
        die "ERROR: $! " if (!$result);

	############################################################

	while(<>)
	{

		if(/^#/)
		{
			print;
			next
		}
		my @F=split;
		next if($F[3]=~/N/);

		$F[8]=~/^GT:DP:AD/ or die "ERROR: $_";
		$F[9]=~/^(.+?):(.+?):\d+,(.+?):/ or die "ERROR: $_";
		my ($GT,$DP,$AF)=($1,$2,int(1000*$3/$2+.5)/1000);

		$F[8]="GT:DP:AF";
		$F[9]="$GT:$DP:$AF";

		if(length($F[3]) eq length($F[4]))
		{
			foreach my $i (0..length($F[3])-1)
			{
				if(substr($F[3],$i,1) ne substr($F[4],$i,1))
				{
					my $key=join "\t",($F[0],$F[1]+$i,$F[2],substr($F[3],$i,1),substr($F[4],$i,1));
					my $val=join "\t",($key,@F[5..@F-1]); $val.="\n";

					if(!$AF{$key} or $AF>$AF{$key})
					{
						$h{$key}=$val;
						$AF{$key}+=$AF;
					}
				}
			}
		}
		else
		{
			my @F3=split //,$F[3];
			my @F4=split //,$F[4];

			while(@F3>1 and @F4>1 and $F3[-1] eq $F4[-1])
			{
				pop @F3;
				pop @F4;
			}
			$F[3]=join "",@F3;
			$F[4]=join "",@F4;

			#######################################################

			my $key=join "\t",@F[0..4];
                        my $val=join "\t",($key,@F[5..@F-1]); $val.="\n";
			if(!$AF{$key} or $AF>$AF{$key})
                        {
				$h{$key}=$val;
                                $AF{$key}+=$AF;
                        }
		}
	}
	my @keys=sort keys %h;
	foreach my $key (@keys)
	{
		$h{$key}=~/(.+):/;
		print "$1:$AF{$key}\n";
	}
}

