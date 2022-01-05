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
	my %h;
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
		if(length($F[3]) eq length($F[4]) and length($F[3])>1)
		{
			foreach my $i (0..length($F[3])-1)
			{
				my $key=join "\t",($F[0],$F[1]+$i,$F[2],substr($F[3],$i,1),substr($F[4],$i,1));
				my $val=join "\t",($key,@F[5..@F-1]); $val.="\n";
				$h{$key}=$val;
			}
		}
		else
		{
			if(length($F[3])>1 and length($F[4])>1)
			{
				my $min=(length($F[3])<length($F[4]))?length($F[3]):length($F[4]);
				$min--;
				$F[1]+=$min;
				$F[3]=substr($F[3],$min);
				$F[4]=substr($F[4],$min);
			}	
			my $key=join "\t",@F[0..4];
                        my $val=join "\t",($key,@F[5..@F-1]); $val.="\n";
                        $h{$key}=$val;
		}
	}
	my @keys=sort keys %h;
	foreach my $key (@keys)
	{
		print $h{$key};
	}
}

