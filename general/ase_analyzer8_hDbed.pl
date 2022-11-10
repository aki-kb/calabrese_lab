#!opt/bin/perl -w
use strict;

## USAGE: perl ase_analyzer8_hDbed.pl <allelic input> all-chr_mm9_10kb-bin.bed <output>
  ## allelic input file is the output from intersect_reads_snps17.pl

## strand of feature does not matter
## UDP method returns reverse strand of reality, take this into account when parsing.

my($in_snps, $in_coords, $outpref)=@ARGV;


my $snphoa=readsnps($in_snps);
my $asehoh=intersect($snphoa, $in_coords);
printer($asehoh, $outpref);

### SUBS ###

#take snp overlapping readand parse into 
# broken up into 10kb windows
sub readsnps {
    my ($in, $species)=@_;
    open (IN, "$in") or die "Cant open infile $in\n";  
    my %coords;

    <IN>;
     while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $seqid=$array[0];
	my $chr=$array[2];
	my $strand=$array[3];
	my $b6start=$array[4];
	my $genmatch=$array[15];

	my $pos=int ($b6start/10000);
	my $info="$b6start\t$strand\t$genmatch";

	@{$coords{$chr}}[$pos].="$info\n";
	#print "$info\n";
    }
    return \%coords;
    close IN;
}
    
#go through coord file (refseqs) and line by line intersect with snp-overlappin
# reads
sub intersect {
    my ($snps, $in)=@_;
    
    open (IN, "$in") or die "Cant open infile $in\n";  

    my $i;
    my %aseinfo;
    #read coord file
    <IN>;
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $gene_name=$array[3];
	my $chr=$array[0];
	my $txstart=$array[1];
	my $txend=$array[2];
	
	#store name of gene as gensym or refseq if no gensym
	my $finalid;
	$finalid="$gene_name\t$chr\t$txstart\t$txend";


	#store name in hash
	if (exists $aseinfo{$finalid}) {
	    $i++;
	    #print "$line\n";
	}

	$aseinfo{$finalid}{'B6'}="0";
	$aseinfo{$finalid}{'CAST'}="0";

	#define bin
	my $refstart=int ($txstart/10000);
	my $refend=int ($txend/10000);
	#define intervals to search
	my $binstart=$refstart-1;
	my $binend=$refend+1;
	
	#define array of potential interacting snps
	my @seqs;

	#if there is data in any of the bins surrounding the gene,
	# store data in @seqs
	for (my $i=$binstart; $i<=$binend; $i++) {
	    if (defined $snps->{$chr}[$i]) {

		my @bininfo;

		#split the info in the bins by\n
		@bininfo=split(/\n/, $snps->{$chr}[$i]);
		push(@seqs, @bininfo);

	    }
	}
	    

	#now search snps within @seqs for overlap
	foreach my $snpline (@seqs) {

	    my @snpinfo=split(/\t/, $snpline);
	    my $pos=$snpinfo[0]; #snp location
	    my $readstrand=$snpinfo[1];
	    my $genmatch=$snpinfo[2];

	 
	    #count all overlapping snps
	    if ($txstart<=$pos && $txend>=$pos) {
		$aseinfo{$finalid}{$genmatch}+=1;
	    }
	    
	}

    }
    #print "$i dups\n";
    return \%aseinfo;
}

sub printer{
    my ($aseh, $outname)=@_;
    open (OUT, ">$outname") or die "cant open outfile name $outname\n";
    
    #print OUT"id\tchr\tstrand\tstart\tend\tb6\tcast\tb6rpm\tcastrpm\n";

    foreach my $coord (keys %{$aseh}) {
	my $b6=$aseh->{$coord}{'B6'};
	my $cast=$aseh->{$coord}{'CAST'};
	
	#my $lb6rpm=$b6/$tr*1000000;
	#my $lcrpm=$cast/$tr*1000000;
	#my $sb6rpm=sprintf("%.2f", $lb6rpm);
	#my $scrpm=sprintf("%.2f", $lcrpm);

	print OUT "$coord\t$b6\t$cast\n";
	

    }
}
