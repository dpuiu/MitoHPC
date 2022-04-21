#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my (%opt,%len);
	#$len{ref}=   3217346917;	#hs38DH size (3355 sequences) 
	#$len{ref}=  3099922541;	#hs38DH size (195 sequences)
	#$len{female}=3160119502;	#hs38DH  based estimates
	#$len{male}=  3110712762;	#hs38DH  based estimates

	$len{ref}=3031865587;		#CHM13 v1.1.	94.23% of prev estimate 
	$len{female}=3054815472;	#CHM13 v1.1     96.66% of prev estimate 
	$len{male}=3008915703;		#CHM13 v1.1;  	96.72% of prev estimate 

	$opt{chrM}=  16569;

	my $result = GetOptions(
		"ref=i"	=>	\$opt{ref},
		"female" =>   \$opt{female},
		"male" =>     \$opt{male},
        );
        die "ERROR: $! " if (!$result);

	if($opt{female})   { $opt{ref}=$len{female} }
	elsif($opt{male})  { $opt{ref}=$len{male}   }
	elsif(!$opt{ref})  { $opt{ref}=$len{ref}    }

	while(<>)
	{
		#0		  1	     2	        3
		#Run              all        mapped     chrM
		#MH0162792.final  740589366  739237125  487382
		#MH0162809.final  763658318  762297733  495743

		my @F=split;
		if(@F>=4)
		{
			if($.==1 and (/^Run/ or /^sample/))
			{
				push @F,"M";
			}
			else
			{
				my $M=($F[3]*2*$opt{ref})/($F[2]*$opt{chrM});
				#$M=int($M*100+.5)/100;
				$M=int($M+.5);
				$F[4]=$M ;
			}
		}
		print join "\t",@F;  
		print "\n";
	}
	exit 0;
}


