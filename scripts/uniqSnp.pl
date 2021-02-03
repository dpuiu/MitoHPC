#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	my %h;
        my $result = GetOptions();
	die "ERROR: $! " if (!$result);

	while(<>)
	{
		next if(/^$/);
		if(/^#/) { print; next }

		my @F=split;
		my $key="$F[1] $1"  if(/.+SM=(.+?);/);
		print unless($h{$key});
		$h{$key}=1;
	}
	exit 0;
}
