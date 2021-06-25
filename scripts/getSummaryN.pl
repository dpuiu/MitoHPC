#!/usr/bin/perl -w
 
use strict;
use warnings;
use Getopt::Long;


###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %options;
	my @keys;
	my %vals;

	# validate input parameters
        my $result = GetOptions(
		"p=s"	=>	\$options{prefix}	
	);
	die "ERROR: $! " if (!$result);
	########################################

        #Run                            H       h       S       s       I       i       Hp      hp      Sp      sp      Ip      ip	A
        #GTEX-1GPI7-0002-SM-DLIOU       17      3       16      0       1       3       2       3       1       0       1       3	20

	while(<>)
	{
		chomp;
		if($.==1)
		{
			@keys=split;
			shift @keys;
		}	
		else
		{
			my @F=split;
			shift @F;
	
			foreach my $i (0..@F-1)
			{
				push @{$vals{$keys[$i]}},$F[$i];
			}		
		}		 
	}
		
	#########################################	
	
	print "$options{prefix}\t" if($options{prefix});
	print join "\t", ("id","count","nonZero","min","max","median","mean","sum");	
	print "\n";
		
	foreach my $key (@keys)
	{	
		my @vals=sort {$a<=>$b} @{$vals{$key}};		
		my $count=scalar(@vals);
		
		my $sum=0;
		my $nonZero=0;
		foreach my $val (@vals) 
		{ 
			$sum+=$val ;
			$nonZero++ if($val);
		}
							
		my $min=$vals[0];
		my $max=$vals[-1];
	
		my $index=int(scalar(@vals)/2);
		my $median=$vals[$index];
                my $mean=int($sum*100/$count+.5)/100;
		
		print "$options{prefix}\t" if($options{prefix});
		print join "\t", ($key,$count,$nonZero,$min,$max,$median,$mean,$sum);
		print "\n"
	}		
	

	exit 0;
}
