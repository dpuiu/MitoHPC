#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	$opt{ref}=   3217346917;	#hs38DH size (3355 sequences) 
	#$opt{ref}=  3099922541;	#hs38DH size (195 sequences)
	$opt{female}=3160119502;
	$opt{male}=  3110712762;
	$opt{chrM}=  16569;

	my $result = GetOptions(
		"ref=i"	=>	\$opt{ref},
		"female=i" =>   \$opt{female},
		"male=i" =>     \$opt{male},
        );
        die "ERROR: $! " if (!$result);

	if($opt{female})  { $opt{ref}=$opt{female} }
	elsif($opt{male}) { $opt{ref}=$opt{male}   }

	while(<>)
	{
		#0		  1	     2	        3	4	
		#Run              all        mapped     chrM    filter
		#MH0162792.final  740589366  739237125  487382  205095
		#MH0162809.final  763658318  762297733  495743  205156

		my @F=split;
		if(/^Run/ or /^sample/)
		{
			push @F,"M";
		}
		else
		{
			my $M=($F[3]*2*$opt{ref})/($F[2]*$opt{chrM});
			#$M=int($M*100+.5)/100;
			$M=int($M+.5);
			#$M="." if($M>10000);
			push @F,$M ;
		}

		print join "\t",@F;  
		print "\n";
	}
	exit 0;
}


