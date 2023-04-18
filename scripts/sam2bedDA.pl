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
	my %len;

	# validate input parameters
	my $result = GetOptions(
		"min_ins=i"	=> \$options{min_ins},
		"max_ins=i"     => \$options{max_ins},
	);
        die "ERROR: $! " if (!$result);

	####################################

	while(<>)
	{
		my @F=split;
                next unless(@F);

		#0					1	2	3	4	5	6	7	8	9
		#L7C2CCXX160303:1:1110:11667:61784	163	chrM	1	26	41S110M	=	59	209	TAA...	*	AS:i:105 XS:i:94 SA:Z:chrM,16529,+,41M110S,26,0;	XA:Z:chr17,+22521368,8S134M9S,8;	MC:Z:151M	MQ:i:60	MD:Z:72A37	NM:i:1	RG:Z:HL7C2.1

		if(/^\@SQ\s+SN:(\S+)\s+LN:(\S+)/)
		{
			$len{$1}=$2;
			next;
		}
		elsif(/^\@/ or @F<11)
                {
                        next;
                }

		next if($F[1]!~/^\d+$/);
		next if($F[1] & 0x4);
		next if($options{mated} and !($F[1] & 0x2) );
		next if($F[2] eq "*");
		next if($F[1] and $F[1] & 0x100) ;  #secondary alignment
		next if($F[8]!~/\d+$/);
		next if($options{min_ins} and abs($F[8])<$options{min_ins});
		next if($options{max_ins} and abs($F[8])>$options{max_ins});

		#################################################

		#SA:Z:chrM,6577,-,101S50M,60,0;
		#if(/\tSA:Z:($F[2],\S+)/)
		if(/\tMC:Z:(\S+)/ )
		{
			my $sacigar=$1;

			#################################################
			my ($ref,$begin,$end,$CIGAR,$cigar,$qry,$strand,$score);

	                $qry=$F[0];
	                #$qry=$qry."/1" if($F[1] & 0x40);
        	        #$qry=$qry."/2" if($F[1] & 0x80);
                	$strand=($F[1] & 0x10 )?"-":"+";
			$ref=$F[2];

	                $CIGAR=$cigar=$F[5];
        	        $begin=$end=$F[3]-1 ;
                	while($cigar and $cigar=~/(\d+)(\w)(.*)/)
                	{
                        	$end+=$1 if($2 eq "M" or $2 eq "D" or $2 eq "N") ;
                        	$cigar=$3;
               	 	}
                	$score=$end-$begin;

			#################################################

			my ($saref,$saqry,$sabegin,$saend,$saCIGAR,$sastrand,$sascore);

			$saqry=$F[0];
                        #$saqry=$saqry."/1" if($F[1] & 0x80);
                        #$saqry=$saqry."/2" if($F[1] & 0x40);
                        $sastrand=($F[1] & 0x20 )?"-":"+";
			if($F[6] eq "=") { $saref=$ref;}else {$saref=$F[6]}
			$sabegin=$saend=$F[7]-1;
			$saCIGAR=$sacigar;

		        while($sacigar and $sacigar=~/(\d+)(\w)(.*)/)
                	{
                        	$saend+=$1 if($2 eq "M" or $2 eq "D" or $2 eq "N") ;
                       		$sacigar=$3;
                	}
			$sascore=$saend-$sabegin;

			if($begin<$sabegin)
			{
				print join "\t",($ref,$begin,$end,$qry,$score,$strand,$CIGAR);print "\t";
				print join "\t",($saref,$sabegin,$saend,$saqry,$sascore,$sastrand,$saCIGAR);print "\n";
			}
			else
                        {
				print join "\t",($saref,$sabegin,$saend,$qry,$sascore,$sastrand,$saCIGAR);print "\t";
                                print join "\t",($ref,$begin,$end,$saqry,$score,$strand,$CIGAR);print "\n";
                         }
		}
	}
	exit 0;
}

