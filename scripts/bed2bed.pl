#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

# help info
my $HELPTEXT = qq~
Program that c...

Usage: $0 files [options]
	  
  INPUT:   
  
  files:	
  
  options:
 
  OUTPUT:  
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %opt;
	$opt{min}=0;
	$opt{add}=0;
	my $result = GetOptions(
		"min=i"	=>	\$opt{min},
		"add=i"	=>	\$opt{add},
		"offset" =>     \$opt{offset}
        );
        die "ERROR: $HELPTEXT " if (!$result);

	while(<>)
	{
		next if(/^$/ or /^#/);

		my @F=split;
		die "ERROR: $_" if(@F<3);

                if($F[0]=~/^(.+):(\d+)-(\d+)$/ && $opt{offset})
                {
                        $F[0]=$1;
                        $F[1]+=$2;
                        $F[2]+=$2;
                }

		$F[1]+=$opt{add};
		$F[2]-=$opt{add};

		push @F,"$F[0]:$F[1]-$F[2]" if(@F==3);
 
		push @F,$F[2]-$F[1] if(@F==4);
		push @F,"." if(@F==5); 

		($F[1],$F[2],$F[5])=($F[2],$F[1],"-") if($F[1]>$F[2]);
		next if(abs($F[4])<$opt{min});
		print join "\t",@F[0..5]; print "\n";
	}
	exit 0;
}
