#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

MAIN:
{
	# define variables
	my %opt;
	my %h;
	my %suspicious;

	my $result = GetOptions(
		"in=s"	=>	\$opt{in},
		"suspicious=s"  => \$opt{suspicious}		
	);
        die "ERROR: $! " if (!$result);

	#######################################################

        if($opt{suspicious})
        {
                open(IN,$opt{suspicious}) or die "ERROR: $!";
                while(<IN>)
                {
                        my @F=split;
                        $suspicious{$F[0]}=1 if(@F);
                }
                close(IN)
        }


        my $AN=0;
	open(IN,$opt{in}) or die "ERROR1: $!";
        while(<IN>)
        {
		chomp;
		next if(/^#/ or /^$/);
		my @F=split;
		next if($suspicious{$F[0]});
		$AN++;
        }
	close(IN);

	#######################################################
	my ($key,$sample,@keys,%keys,%GT_DP_AF);
	while(<>)
	{
		#chrM	3243	rs199474657	A	G	.	PASS	SM=1121939_23193_0_0;TRN=TRNL1;CR=0.740946;CP=9.093	GT:DP:AF	0/1:993:0.146

		#chr   pos  ID	ref  alt     filters  AC_hom  AC_het  AF_hom        AF_het        AN     max_observed_heteroplasmy
		#chrM  3    .	T    C       PASS     19      1       3.3667646E-4  1.7719814E-5  56434  9.9700e-01

		if(/^#/)
		{
		}
		else
		{
			my @F=split;
			my $key="$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]";

			if($F[7]=~/^SM=.+?;(.+)/)
			{
				$h{$key}{filter}=$1;
			}
			else
			{
				$h{$key}{filter}="."
			}

			if(/.+:1$/)
			{
				$h{$key}{AC_hom}++;
			}
			elsif(/.+:(0\.\d+)$/)
                        {
                                $h{$key}{AC_het}++;
				my $AF=$1;
				$h{$key}{max_observed_heteroplasmy}=$AF if(!$h{$key}{max_observed_heteroplasmy} or $AF>$h{$key}{max_observed_heteroplasmy});
			}
			elsif(/.+:0$/)
			{
				next
			}
			else
			{
				die "ERROR: $_"
			}

  			if(!$keys{$key})
                       	{
                               	$keys{$key}=1;
                               	push @keys,$key;
                        }

		}
	}

	print "chr\tpos\tID\tref\talt\tfilters\tAC_hom\tAC_het\tAF_hom\tAF_het\tAN\tmax_observed_heteroplasmy\n";
	foreach my $key (@keys)
	{
		my @F;

		foreach ("AC_hom","AC_het","max_observed_heteroplasmy")
		{
			$h{$key}{$_}=0 unless(defined($h{$key}{$_}));
		}

                $h{$key}{AF_hom}=int(200000*$h{$key}{AC_hom}/$AN+.5)/200000;
		$h{$key}{AF_het}=int(200000*$h{$key}{AC_het}/$AN+.5)/200000;


		push @F,($key,$h{$key}{filter},$h{$key}{AC_hom},$h{$key}{AC_het},$h{$key}{AF_hom},$h{$key}{AF_het},$AN,$h{$key}{max_observed_heteroplasmy});
		print join "\t",@F; 
		print "\n";
	}
	exit 0;
}

