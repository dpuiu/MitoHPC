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
		elsif($.==1 or /^#/)
		{
			print;
			next;
		}

		chomp;
		my @F=split /\t/;

		my $F3=substr($chrM,$F[1]-1,length($F[3]));
		if($F[3] ne $F3)
		{
			if(length($F[3]) eq length($F[4])) { next  }
			else                               { $F[3]=$F3 }
		}	

		print join "\t",@F; print "\n";	
	}
}

