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
		"ref_len=s"	=>	\$options{ref_len}
	);
	die "ERROR: $! " if (!$result);
	die "ERROR: missing reference .fai"    unless($options{ref_len});

	open(IN,$options{ref_len}) or die;
        while(<IN>)
	{
		chomp;
                next if(/^#/ or /^$/);
             	next unless(/^(\S+)\s+(\d+)/);

        	$len{$1}=$2;
	}
	close(IN);
	########################################
	while(<>)
	{
		my @F=split /\t/;
		if(/^\@SQ\s+SN:(\S+)\s+LN:(\S+)/)
		{
			if($len{$1})
			{
				$F[2]="LN:$len{$1}\n";
			}
		}
		elsif(/^\@/)
                {
                }
		else
		{
			if($F[7])
			{
				if($F[6] eq "=" and $len{$F[2]} and $F[7]>$len{$F[2]})
				{
					$F[7]%=$len{$F[2]}
				}
				elsif($len{$F[6]} and $F[7]>$len{$F[6]})
                               	{
                                       	$F[7]%=$len{$F[6]}
                               	}
			}

			if($F[3] and $len{$F[2]} and $F[3]<=$len{$F[2]})
			{
				my ($cigar1,$cigar2)=("","");
				my ($trim1,$trim2)=(0,0);
				my $pos=$F[3]-1;

				if($F[5]=~/^(\d+)([SH])(.+)/) 
				{
					$cigar1=$1.$2;
					$trim2+=$1;
					$F[5]=$3;
				}

				while($F[5]=~/^(\d+)(\w)(.*)/)
				{
					if($pos+$1<=$len{$F[2]})
					{
						$cigar1.=$1.$2;
						$trim2+=$1 unless($2 eq "D");
					}
					elsif(!$trim1 and $2 eq "D")
					{
						$cigar1.=$1.$2;
					}
					elsif(!$trim1)
					{
						$trim1=$pos+$1-$len{$F[2]};
						my $keep1=$1-$trim1;

						if($keep1)
						{
							$trim2+=$keep1;
							$cigar1.=$keep1.$2;
						}

						$cigar2=$trim2."S".$trim1.$2;
					}
					else
					{
						$cigar2.=$1.$2;
						$trim1+=$1 unless($2 eq "D");
					}

					$pos+=$1 unless($2 eq "D");
					$F[5]=$3;
				}

				if($trim1)
				{
					$cigar1.=$trim1."S";
				}

				if($cigar1=~/^(.+?)(\d+)S(\d+)S$/)
				{
					$cigar1=$1.($2+$3)."S";
				}

				if($cigar2=~/M/) 
				{ 

					print join "\t",($F[0],$F[1]+2048,$F[2],1,$F[4],$cigar2,@F[6..scalar(@F)-1]);
				}

			        if($cigar1=~/M/) 
                               	{ 
                                  	$F[5]=$cigar1;
                               	}
				else
				{
					@F=()
				}

			}
			elsif($F[3] and $len{$F[2]} and $F[3]>$len{$F[2]})
			{
				$F[3]=$F[3]%$len{$F[2]};
			}
		}
		print join "\t",@F if(@F); 	
	}
	exit 0;
}

