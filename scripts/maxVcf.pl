#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that computes each SNP major allele 

        EXAMPLE:
                 cat I.vcf | maxVcf.pl
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my (%max,%line,%sum,%af,%line2,%af2,%line3);
	while(<>)
	{
		if(/^#/) { print; next; }

		chomp;
		my @F=split /\t/;

		my $POS=$F[1];

		my $AF=1;
		$AF=$1 if($F[9]=~/.+:(\S+)/);
                if(!$max{$POS} or $AF>$max{$POS})
                {
                        $max{$POS}=$AF;
                        $line{$POS}=$_;
			$af{$POS}=$AF;
                }

		$sum{$POS}+=$AF;
	}


	foreach my $POS ( sort {$a <=> $b} keys %line )
	{
		if($max{$POS}>1-$sum{$POS})
		{
			$line2{$POS}=$line{$POS};
			$af2{$POS}=$af{$POS};
		}
	}


        foreach my $POS ( sort {$af2{$b} <=> $af2{$a}} keys %line2 )
	{
		#print "\n#$line2{$POS}\n";

		my @POS= sort {$a <=> $b} keys %line3;
		while(@POS and $POS > $POS[0])
		{
			shift @POS
		}

		if(@POS)
		{
			#print "#$line3{$POS[0]}\n";
			my @P=split /\t/,$line2{$POS};

			if($POS[0]-$POS>length($P[3])-length($P[4]))
			{
				$line3{$POS}=$line2{$POS};
                                #print "#ok\n";
			}
		}
		else
		{
			$line3{$POS}=$line2{$POS};
			#print "#ok\n";
		}
	}

	foreach my $POS ( sort {$a <=> $b} keys %line3 )
	{
		print $line3{$POS},"\n"
	}
        
}

