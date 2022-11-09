#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Time::Local;

## USAGE: perl bigbowtie_to_wig3_mm9.pl <input>.sam <outfile> <color>

###   MAIN   ###

my ($in_bow, $name, $color)=@ARGV;

my $headerName = $name;
$headerName =~ s/\A.*\///g;
my $outName = $name.".wig";

my ($second,$minute,$houre,@reste) =localtime(time);
print $houre.":".$minute.":".$second." start making wigs...\n";

my %colors = (
	'pink' => '255,181,197',
	'babyblue' => '135,206,250',
	'grey' => '136,136,136',
	'black' =>'0,0,0',				
	'red' =>'255,0,0',				
	'green' =>'0,51,0',	
	'blue' =>'0,0,255',				
	'yellow' =>'255,204,0',
	'lightgreen' =>'0,112,0',
	'purple' =>'51,0,51',			
	'navyblue'=>'0,0,51',
	'maroon'=> '51,0,0',
	'brown'=>'51,26,0',				
	'orange'=>'255,77,0',
	'magenta'=>'255,0,255',
	);
unless (exists $colors{$color}) 
{
	die "$color is not an accepted color.\nPlease use pink, babyblue, grey, black, red, green, blue, yellow, lightgreen, purple, navyblue maroon, brown, orange, or magenta.\n";
}																															

my %lengths = (					## based on mm9
	'chr1' => '197195432',
	'chr2' => '181748087',
	'chr3' => '159599783',
	'chr4' => '155630120',
	'chr5' => '152537259',
	'chr6' => '149517037',
	'chr7' => '152524553',
	'chr8' => '131738871',
	'chr9' => '124076172',
	'chr10' => '129993255',
	'chr11' => '121843856',
	'chr12' => '121257530',
	'chr13' => '120284312',
	'chr14' => '125194864',
	'chr15' => '103494974',
	'chr16' => '98319150',
	'chr17' => '95272651',
	'chr18' => '90772031',
	'chr19' => '61342430',
	'chrX' => '166650296',
	#'chrY' => '15902555',
	);

open (OUT, ">".$outName) or die "Cannot open ".$outName."!\n";
print OUT "track type=wiggle_0 visibility=full" .
	" name=\"".$headerName."\"  color=".$colors{$color}.
	" maxHeightPixels=128:40:11 ".
	" group=\"user\" graphType=bar".
	" priority=100 viewLimits=0:2.8 autoscale=on\n";

foreach my $chr (sort keys %lengths) {

	my $chromolen = $lengths{$chr};
	my $wighoa = &readbow($in_bow, $chr, $chromolen);

	if ($wighoa !~ /nochr/) {
		my ($sec,$min,$hour,@rest) = localtime(time);
		print $hour.":".$min.":".$sec." print wig for ".$chr."...\n";
		printer($wighoa);
	}
}

### SUBS ###
## store bed hits in an array that counts 1 for every 50bp covered by 
## hit.

sub readbow {
	my ($in, $refchr, $clen)=@_;
	my %coords;

	my $togrep = $refchr."[[:space:]]";
	my $temp = $name."_".$refchr."_temp.txt";
	
	print "grep ".$togrep." ".$in." ".$temp."\n";
	`grep $togrep $in > $temp`;
	#my $wc=`wc temp.txt`;
	#print "$wc lines in temp.txt\n";

	open (IN, $temp) or die "Cannot open the temporary file ".$temp."!\n";

	my $j;
	my $k;
	while (my $line=<IN>) {
		chomp $line;
		my @array=split(/\t/, $line);

		my $chr=$array[2];
		my $start=$array[3];
		my $end=$start+300;

		if ($end>=$clen) {
			$end=$clen-1;
		}

		if ($refchr !~ /^$chr$/) {
			print $refchr." does not match:\n".$line."\n";
		}

		if ($start !~ /[a-z]/) {

			## bins to add seqs to
			my $fbin = int ($start/50);
			my $ebin = int ($end/50);
		
			## counter to add seq to bins
			my $i = $fbin;
		
			## add a count to each 50mer bin of the chr@ between start and
			## end
			while ($i<=$ebin) {
				@{$coords{$chr}}[$i]+=1;
				$i++;
				$k++;
			}
		}
		$j++;
	}

	close IN;
	`rm -f $temp`;

	if (defined $j) {
		print $j." lines read, ".$k." total hits from ".$refchr."\n";
		return \%coords;
	} else {
		return 'nochr';
	}
}

## sub to print out wiggle format
sub printer {
	my ($counts)=@_;

	foreach my $chr (sort keys %{$counts}) {
		print OUT "variableStep chrom=".$chr." span=50\n";
		my @hits = @{$counts->{$chr}};
		my $i;
		foreach my $hit (@hits) {
			if (defined $hit) {

				my $longlog=(log $hit)/(log 10);
				my $log=sprintf("%.2f", $longlog);
				my $pos=50*($i);
				#print "0$hit $longlog $log\n";
				print OUT $pos." ".$hit."\n";
			}
			$i++;
		}
	}
}
