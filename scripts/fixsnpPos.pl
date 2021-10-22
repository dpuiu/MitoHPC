#!/usr/bin/env perl
 
use strict;
use warnings;
use Getopt::Long;

my $HELP = qq~
Program that recalculates snp positions

        EXAMPLE:

		....
~;

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my %opt;
	my $result = GetOptions(
		"ref=s"		=> \$opt{ref},
		"rfile=s"       => \$opt{rfile},
 	        "file=s" 	=> \$opt{file},
 	);
	die "ERROR: $! " if (!$result);
	die "ERROR: $! " if (!defined($opt{ref}) or !defined($opt{rfile}) or !defined($opt{file}));

	#####################################################################################

	my %diff;
	open(IN,$opt{file}) or die "ERROR:$!";
        while(<IN>)
        {
		if(!/^#/ and /INDEL/)
		{
			my @F=split;
			$diff{$F[1]}=length($F[4])-length($F[3]);
		}
        }
	close(IN);

	my $MT="";
	open(IN,$opt{rfile}) or die "ERROR:$!";
        while(<IN>)
        {
                if(/^>/){}
		else
		{
			chomp;
			$MT.=$_;
		}
	}

	if($opt{ref} eq "rCRS" or $opt{ref} eq "chrM") 
	{
		$diff{3107}=-1
	}

	####################################################################################
	my %diff2;
	my $diff=0;
	my @pos=(sort {$a<=>$b} keys %diff);
	foreach my $pos (@pos)
	{
		$diff+=$diff{$pos};
		#$diff{$pos}=$diff; 		# Dec 8
		$diff2{$pos+$diff}=$diff	# Dec 8
	}
	my @pos2=(sort {$a<=>$b} keys %diff2);

	###################################################
	##reference=file://out/chrM.L0/chrM.L0.mutect2.fa>
	##contig=<ID=chrM.L0,length=16567>

	while(<>)
	{
		if(/^$/)
		{
			next;
		}
		elsif(/^##reference/ or /^##contig/)
		{
			print "##reference=file://$opt{rfile}\n";
			print "##contig=<ID=$opt{ref},length=16569>\n" if($opt{ref} eq "rCRS" or $opt{ref} eq "chrM");
		}
		elsif(/^#/)
		{
			print;
			next;
		}
		else
		{
			chomp;
			my @F=split /\t/;
			my $diff2=0;
			my $i=0;
			while($i<@pos2 and $F[1]>$pos2[$i])
			{
				$diff2=$diff2{$pos2[$i]};
				$i++
			}

			$F[0]=$opt{ref};
			$F[1]-=$diff2; 
			#$F[1]+=$diff2;

			my $F3=substr($MT,$F[1]-1,length($F[3]));	# new 02/24/2021
			if($F3 eq $F[4] and $F[-2]=~/^GT:AD:AF/ and $F[-1]=~/([01])(.)([01]):(\d+),(\d+):(0\.\d+)(:.+)$/)
			{
				$F[4]=$F[3];
				#$F[-1]=(1-$1).$2.(1-$3).":$5,$4:".(1-$6).$7."***";

				my $F_1="";
				while($F[-1]=~/(.+?):(.+)/)
				{
					$F[-1]=$2;

					if($1=~/^(0\.\d+)$/ or $1=~/^(\d+.*e-\d+)$/)      { $F_1.=(1-$1).":" }
					elsif($1=~/(\d+)(\D)(\d+)/)    { $F_1.="$3$2$1:"  }
					else                           { $F_1.="$1:"      }
				}
				$F_1.=$F[-1];
				$F[-1]=$F_1;
			}
			$F[3]=$F3;

			print join "\t",@F;
			print "\n";
		}
	}
}

