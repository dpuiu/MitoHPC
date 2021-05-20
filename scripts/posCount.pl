#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;

	my $result = GetOptions(
	"total=i"	=>	\$opt{total}
	);
        die "ERROR: $! " if (!$result);

	my %pos;
	while(<>)
	{
		next if(/^#/);
		my @F=split;
		my $t=(/AF=0/ or /:0\.\d+$/)?"h":"H";

		$pos{"$F[0]\t$F[1]"}{$t}++;
		$t.=(/INDEL/)?"I":"S";

		$pos{"$F[0]\t$F[1]"}{$t}++;
	}

	print join "\t",("#chr","pos","H","h","HS","hS","HI","hI"); print "\n";
	foreach my $k (keys %pos)
	{
		foreach my $t ("H","h","HS","hS","HI","hI")
		{
			$pos{$k}{$t}=0 unless($pos{$k}{$t});
		}

                foreach my $t ("H","h","HS","hS","HI","hI")
                {
			
                        $pos{$k}{$t}=int(10000*$pos{$k}{$t}/$opt{total}+.5)/100 if($opt{total});
                }

		print join "\t",($k,$pos{$k}{"H"},$pos{$k}{"h"}, $pos{$k}{"HS"},$pos{$k}{"hS"},$pos{$k}{"HI"},$pos{$k}{"hI"});
		print "\n";
	}
	exit 0;
}

