#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that corrects mutserve output

        EXAMPLE:
		cat I.mutserve.vcf | fixmutserveVcf.pl -file rCRS.fa
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %opt;
	my $result = GetOptions(
 	               "file=s" 	=> \$opt{file},
 	       );
	die "ERROR: $! " if (!$result);
	die "ERROR: $! " if (!defined($opt{file}));

	my $chrM="";
	open(IN,$opt{file}) or die "ERROR:$!";
        while(<IN>)
        {
                if(/>/) {}
                else
                {
                        chomp;
                        $chrM.=$_;
                }
        }
	close(IN);
	###################################################

	my @P;
	while(<>)
	{

		if(/^$/)
		{
			next;
		}
		elsif(/^#/)
		{
			print;
			next;
		}

		chomp;
		my @F=split /\t/;

		if($F[1]==3105    and $F[4] eq "*") { next }
		elsif($F[1]==3106 and $F[4] eq "*") { next }
		elsif($F[1]==3107)                  { next }
		elsif($F[1]==3108 and $F[4] eq "*") { next }
		elsif(@P and $P[4] eq "*")
		{
			if($P[0] eq $F[0] and $F[1]-$P[1]==1 and $F[4] eq "*")
	               	{
        	               	$P[3].=$F[3]
               		}
			else
			{
				$P[1]--;
				$P[4]=substr($chrM,$P[1]-1,1);
				$P[3]="$P[4]$P[3]";
				print join "\t",@P;print "\n";
				@P=();
			}
		}
		elsif($F[4] eq "*")
		{
			@P=@F
		}
		else
		{
			print;print "\n";
			@P=();
		}
	}

	if(@P)
	{
		$P[1]--;
		$P[4]=substr($chrM,$P[1]-1,1);
		$P[3]="$P[4]$P[3]";

		print join "\t",@P; 
		print "\n"
	}
}

