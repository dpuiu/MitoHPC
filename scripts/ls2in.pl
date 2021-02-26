#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	$opt{out}="out";

	my $result = GetOptions(
		"out=s"	=>	\$opt{out}
	);
        die "ERROR: $! " if (!$result);

	while(<>)
	{
		chomp;
		next unless(/\.bam$/ or /\.cram$/);
		my @F=split;

		print "$1\t$F[-1]\t$opt{out}/$1/$1\n" if($F[-1]=~/.+\/(\S+)\./ or $F[-1]=~/(\S+)\./);
		#print "$1\t$F[-1]\t$opt{out}/$1/$1\n" if($F[-1]=~/.+\/(\S+?)\./ or $F[-1]=~/(\S+?)\./);
	}
	exit 0;
}

