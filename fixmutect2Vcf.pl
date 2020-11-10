#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
        my %opt;
        my $result = GetOptions(
                       "file=s"         => \$opt{file},
               );
        die "ERROR: $! " if (!$result);

	############################################################

	my @P;
	my $PAF;
	while(<>)
	{
		if(/^#/) 
		{ 
			print; next 
		} 

		my @F=split;
		if(@P and $P[1]==$F[1] and $P[4] eq $F[4] and $P[7]=~/AF=/ and $F[7]=~/AF=/ ) 
		{ 
			$P[7]=~/AF=(\S+)/; $PAF=$1; $F[7]=~/(.*AF=)(\S+)/; $F[7]=$1.($PAF+$2) 
		} 
		elsif(@P)
		{ 
			print join "\t",@P; print "\n"
		} 

		@P=@F
	}

	if(@P) 
	{ 
		print join "\t",@P; 
		print "\n"
	}
}

