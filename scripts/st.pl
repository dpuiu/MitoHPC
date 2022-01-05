#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;


    MAIN:
    {
        # define variables
        my %opt;
	my (@c,$s,$n);

        my $result = GetOptions();
        die "ERROR: $! " if (!$result);

        while(<>)
        {
		chomp;
		/^\d+$/ or /^\d+\.\d+$/ or die "ERROR $_";
		push @c,$_;
		$s+=$_;
        }

	@c=sort {$a<=>$b} @c;

	if(@c)
	{
		$n=@c;
		print join "\t",("count","min","q1","median","q3","max","mean"); print "\n";
		print join "\t",(scalar(@c),$c[0],$c[int($n/4+.5)],$c[int($n/2+.5)],$c[int($n*3/4+.5)],$c[-1],int($s*100/$n+.5)/100); print "\n";
	}

	exit 0;
    }

