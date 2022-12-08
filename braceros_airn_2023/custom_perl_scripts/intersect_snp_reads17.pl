#!opt/bin/perl -w
use strict;
use Getopt::Std;

## USAGE: perl intersect_reads_snps17.pl <b6 input>.sam <cast input>.sam sanger_mm9_09_14_11 <n/y, paired end?> <output>

my($in_B6sam, $in_CASTsam, $in_snps,  $peflag, $outpref)=@ARGV;

## pe flag must be y or n
if (($peflag ne "y") && ($peflag ne "n")) {
    die "pe flag must be y or n\n";
}

my %species_hash=("B6"   =>"$in_B6sam", 
		  "CAST" =>"$in_CASTsam" );

my $snphoa=readsnps($in_snps);


## for each genome alignment file, find snp overlapping reads
foreach my $species_name (keys %species_hash) {
    print "$species_name\n";
    intersect($snphoa, $species_hash{$species_name}, $outpref, 
	    $species_name, $peflag);
}

my $finalname=parsedups($outpref);
summary($finalname);

### SUBS ###

## take sanger informative snp/indel file and file according to chromosomes 
## broken up into 20kb windows
## sam is 1-based, coord info in my file is 0 based
sub readsnps {
    my ($in)=@_;
    open (IN, "$in") or die "Cant open infile $in\n";  
    my %coords;

    <IN>;
 
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $snpid=$array[0];
	my $chr=$array[1];
	my $B6start=$array[2];
	my $B6end=$array[3];
	my $B6seq=$array[4];
	my $caststart=$array[5];
	my $castend=$array[6];
	my $castseq=$array[7];

	## store \n delimited list, by setting the position to
	## $end instead of $start, you account for the 1-based coordinate
	my $pos=int ($B6end/5000);
	
	my $sinfo="$B6end\t$B6seq\t$snpid";
	@{$coords{$chr}{"B6"}}[$pos].="$sinfo\n";
	  
	my $cinfo="$castend\t$castseq\t$snpid";
	@{$coords{$chr}{"CAST"}}[$pos].="$cinfo\n";

    }
    close IN;
    print "snps read\n";
    return \%coords;
}

## go through bowtie file and line by line intersect with snps
sub intersect {
    my ($snps, $in, $out, $species, $pe)=@_;
    
    open (IN, "$in") or die "Cant open infile $in\n";  

    my $outname="$out"."_$species";
    open (OUT, ">$outname") or die "cant open outfile name $outname\n";

    ## file for reads with unread cigar strings
    my $cig="$outname"."_CIG";
    open (CIG, ">$cig") or die "cant open $out cig\n";
    
    my $rc="0"; #reads counted (both mates match)
    my $tr="0"; # total reads in file
    my $sct="0";# reads that were soft clipped by star
    
    ## read sam file
    ## skip headers if they exist

    ## strings to store forward and reverse mate data
    my ($fdata, $rdata, $linecount);

    while (my $fileline=<IN>) {
	chomp $fileline;

	## skip headders
	if (($fileline=~/^\@HD/) ||($fileline=~/^\@SQ/)||($fileline=~/^\@PG/)||($fileline=~/^\@CO/)) {
	    #print "$fileline\n";
	} else {
	    my ($subreads, $pass, $tot, $scp, $id, $chr, $strand, $start,
		$ocigar, $mate, $readspan, $fullseq)=splitsam($fileline, $out, $pe);
	    
	    $rc=$rc+$pass;
	    $tr=$tr+$tot;
	    $sct=$sct+$scp;
	    
	    ## cycle through each fragment of the read and accumulate overlapping snps
	    ## , then want to store read id, chr, strand, f start,
	    ## f end r start r end forward cigar, reverse cigar, genome, f snps, r snps
	    
	    #list of snps that overlap
	    my $snpovers;
	    
	    foreach my $sr (@{$subreads}) {
		my @info=split(/\t/, $sr);
		my $rid=$info[0];
		my $rchr=$info[1];
		my $rstart=$info[2];
		my $rend=$info[3];
		my $rstrand=$info[4];
		my $rseq=$info[5];	
		my $rq=$info[6];
		my $frag=$info[7];
		my $rcigar=$info[8];

		$rchr=~s/.fa//;
		
		#print "$sr\n";
		
		my $frl=length($rseq);
		
		## define bin
		my $ref=int ($rstart/5000);
		## define intervals to search
		my $a1pos=$ref-1;
		my $a2pos=$ref;
		my $a3pos=$ref+1;
		
		## define array of potential interacting snps
		my @seqs;
		
		## if there is data in any of the 3 bins surrounding the gene,
		## count the number of reads that surround the bin
		if (defined $snps->{$rchr}{$species}[$a1pos]  || 
		    defined $snps->{$rchr}{$species}[$a2pos]  || 
		    defined $snps->{$rchr}{$species}[$a3pos]) {
		    
		    my (@a1, @a2, @a3);
		    
		    if (defined $snps->{$rchr}{$species}[$a1pos]) {
			@a1=split(/\n/, $snps->{$rchr}{$species}[$a1pos]);
			push(@seqs, @a1);
		    }
		    
		    if (defined $snps->{$rchr}{$species}[$a2pos]) {
			@a2=split(/\n/, $snps->{$rchr}{$species}[$a2pos]);
			push(@seqs, @a2);
		    }
		    
		    if (defined $snps->{$rchr}{$species}[$a3pos]) {
			@a3=split(/\n/, $snps->{$rchr}{$species}[$a3pos]);
			push(@seqs, @a3);
		}
		    
		    ## now search snps within @seqs for overlap
		    ## stop after the first snp that overlaps.
		    my $flag;
		    foreach my $snpline (@seqs) {
			
			## check to see if there is already an overlapping snp
			if (defined $flag) {
			    last;
			}
			
			my @snpinfo=split(/\t/, $snpline);
			my $pos=$snpinfo[0]; #snp location
			my $snpseq=$snpinfo[1]; #reference seq at loc
			my $snpid=$snpinfo[2]; #line number of snp on original file
			
			## match first overlapping snp
			if ($rstart<=$pos && $rend>=$pos) {
			    
			    my $snppos=$pos-$rstart;
			    
			    ## find seq of overlapping base to ensure it is 
			    ## not a mismatch
			    my @bases=split ('', $rseq);
			    my $rsnpseq=$bases[$snppos];
			    
			    
			    ## if snp matches read seq at that position, we are good
			    if ($rsnpseq=~/$snpseq/) {
				
				## find qual of overlapping base
				my @quals=split ('', $rq);
				my $snpqual=$quals[$snppos];
				my $qscore=ord($snpqual)-33;
				
				## define flag to end foreach loop
				$flag=1;
				
				## only count if qscore is above 20
				if ($qscore>=20) {
				    $snpovers.="$snpid|$snpseq,";
				}
			    } else {
				#print "no match at snp position\n";
			    }
			    
			}
		    }
		}	    
	    }
	    
	    ## now determine if read overlapped snp at all
	    ## my ($subreads, $pass, $tot, $scp, $id, $chr,
	    # $strand, $start, $ocigar, $mate, $readspan)=splitsam($fileline, $out, $pe);
	    if ($pass>0) {
		
		my $end=$start+$readspan-1;
		
		if (!defined $snpovers) {
		    $snpovers="NA";
		} else {
		    chop $snpovers;
		}
		
		#print "$id\t$snpovers\t$mate\t$chr\n";
		
		if ((defined $fdata) && (defined $rdata)) {
		    
		    ## now both are defined, so print the previous to outfile if there is a snp and reset the data
		    my $together="$fdata\t$rdata\t$species";
		    if ($together=~/snpid/) {
			print OUT "$together\n";
		    }
		    
		    #reset and redefine
		    $fdata=undef();
		    $rdata=undef();
		    if ($mate=~/f/) {
			$fdata="$id\t$mate\t$chr\t$strand\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    } elsif ($mate=~/r/) {
			$rdata="$mate\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    } else {
			$fdata="$id\t$mate\t$chr\t$strand\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
			$rdata="$mate\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    }
		    
		} else {
		    
		    if ($mate=~/f/) {
			$fdata="$id\t$mate\t$chr\t$strand\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    } elsif ($mate=~/r/) {
			$rdata="$mate\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    } else {
			$fdata="$id\t$mate\t$chr\t$strand\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
			$rdata="$mate\t$start\t$end\t$ocigar\t$snpovers\t$fullseq";
		    }
		    
		}
	    }
	}
    }

    #final line
    my $together="$fdata\t$rdata\t$species";
    if ($together=~/snpid/) {
	print "$together\n";
    }

    close OUT;
    close IN;
    print "$rc reads considered from $tr reads in $in, $sct reads soft clipped\n";
}

## sub to go through each output file and make sure that seq id's are unique.
## any id that matches both B6 and CAST is ambiguous and not allele specific.
sub parsedups {
    my ($in)=@_;
    
    my (%all, %dups);
    
    ## create a hash containing ids that are in both CAST and B6 files
    ## this also gets marks reads that have snps in more than one location in read
    
    my @species=("B6", "CAST");

    ## go through 129 first
    my $inname="$in" . "_B6";
    open (IN, "$inname") or die "no such file as $inname\n";
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $id=$array[0];	
	$all{$id}=1;
	
    }
    close IN;
    
    ## now go through cast
    $inname="$in" . "_CAST";
    open (IN, "$inname") or die "no such file as $inname\n";
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $id=$array[0];

	if (defined $all{$id}) {
	    $dups{$id}=1;
	}
    }
    close IN;

    
    ## now print final merged file that does not contain dups.
    ## also print one line per read that overlaps a snp
    
    my $outname="$in"."_final";
    open (OUT, ">$outname") or die "cant open final outfile\n";
    print OUT "readid\tmate\tchr\tstrand\tf_start\tf_end\tf_cig\tf_snps\tf_seq".
	"\tmate\tr_start\tr_end\tr_cig\tr_snps\tr_seq\tgen\n";

    #now go through both files again
    foreach my $spe (@species) {
	my $inname="$in" . "_$spe";
	open (IN, "$inname") or die "no such file as $inname\n";
	
	while (my $line=<IN>) {
	    chomp $line;
	    my @array=split(/\t/, $line);
	    my $id=$array[0];

	    if ((!defined $dups{$id})) {
		print OUT "$line\n";
	    }
	}
    }
    close IN;
    close OUT;
    return \$outname;
}

## sub to get summary info from _final file
sub summary  {
    my ($ref)=@_;
    my $in=$$ref;
    open (IN, "$in") or die "no such file as $in\n";

    my ($b6, $cast, $chrXb6, $chrXcast, $total);
    <IN>;
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $chr=$array[2];
	my $match=$array[15];

	$total++;
	if ($match=~/B6/) {
	    $b6++;
	} elsif ($match=~/CAST/) {
	    $cast++;
	}
	
	if ($chr=~/chrX/) {
	    if ($match=~/B6/) {
		$chrXb6++;
	    } elsif ($match=~/CAST/) {
		$chrXcast++;
	    }
	}
    }
    close IN;
    print "total\t$total\nb6\t$b6\nCAST\t$cast\nchrXb6\t$chrXb6\nchrXCAST\t$chrXcast\n";
}

## split each read into multiple read equivalents based on cigar string
sub splitsam {
    my ($line, $out, $pe)=@_;

    my $cig="$out"."_CIG";
    open (CIG, ">$cig") or die "cant open $out cig\n";

    my $i="0";
    my $j="0";
    my $sc="0"; ## counting reads w/soft clips
    
    chomp $line;
    my @array=split(/\t/, $line);
    my $id=$array[0];
    my $flag=$array[1];  
    my $chr=$array[2];
    $chr=~s/.fa//;
    
    my $start=$array[3];  
    my $mq=$array[4];  
    my $cigar=$array[5]; 
    my $seq=$array[9];
    my $qual=$array[10];

    my $ocigar=$cigar;
    my $oseq=$seq;
  
    ## info for each subread
    my @subreads;

    ## forward or reverse read or not paired end
    my $mate="f";
    if ($pe=~/n/) {
	$mate="NA";
    }

    $j++;

    ## total genomic span of read
    my $readspan;

    ## strand
    my $strand;
    
    ## not yet sure what to do with cigars with these operators
	## if it does not match these
    if ($cigar!~/(P|H|X)/) {
	
	## find strand
	## flag webpage:
        #https://broadinstitute.github.io/picard/explain-flags.html

	if ($pe=~/y/) {

	    if ($flag eq 99) { 
		$strand='+';
	    } elsif ($flag eq 147) {
		$strand='+';
		$mate="r";
	    } elsif ($flag eq 83) {
		$strand='-';
	    } elsif ($flag eq 163) {
		$strand='-';
		$mate="r";
	    }
	} else {
	    if ($flag & 16) { 
		$strand='-';
	    } else {
		$strand='+';
	    }
	}
	
	## parsing cigar string

	if (defined $strand) {
	    $i++;
	    ## string for cigar length
	    my $cl=length($cigar);
	    ## array info for cigar
	    my @cinfo;
	    while ($cl>0) {
		$cigar=~s/^(\d+\D)(.*)/$2/;
		push (@cinfo, $1);
		$cl=length($cigar)
	    }
	    
	    ## go through cigar info and make new start end coords, seq for the read
	    my ($end);
	    my $sr;
	    foreach my $action (@cinfo) {
		
		my $srinfo;
		
		if ($action=~/(\d+)M/) {

		    ## M is match -- store this information as sequence
		    $sr++;
		    $end=$start+$1-1;
		    $readspan+=$1;
		    my $nid="$id"."__"."$sr"."$mate";
		    
		    ## extract from beginning of seq
		    my $toextract=$end-$start+1;
		    $seq=~s/^([a-zA-Z]{$toextract})(.*)/$2/;
		    my $sseq=$1;
		    
		    $qual=~s/^([\S+]{$toextract})(.*)/$2/;
		    my $squal=$1;
		    
		    $srinfo="$nid\t$chr\t$start\t$end\t$strand\t$sseq\t$squal\t".
			"$action\t$ocigar";
		    
		    push (@subreads, $srinfo);
		    $start=$end+1;
		    
		} elsif  ($action=~/(\d+)N/) {
		    $start=$start+$1;
		    $readspan+=$1;
		    
		} elsif  ($action=~/(\d+)D/) {
		    
		    ## read missing \d base in reference
		    $start=$start+$1;
		    $readspan-=$1;

		} elsif  ($action=~/(\d+)S/) {
		    
		    ## soft clip -- here the sequence and qual scores need to be clipped but not the
		    ## starting alignment position
		    ## soft clips at the end of the read do not come with any more sequence either.
		    $start=$start;
		    $seq=~s/^([a-zA-Z]{$1})(.*)/$2/;
		    $qual=~s/^([\S+]{$1})(.*)/$2/;
		    $sc++;
		    
		} elsif  ($action=~/(\d+)I/) {
		    
		    ## read has base inserted relative to reference
		    $seq=~s/^([a-zA-Z]{$1})(.*)/$2/;
		    $qual=~s/^([\S+]{$1})(.*)/$2/;
		} 
	    }
	} else {
	    ## print out to file with PHX in cigar
	    print CIG "$line\n";
	}
    }
    return (\@subreads, $i, $j, $sc, $id, $chr, $strand, $start, $ocigar, $mate, $readspan, $oseq);
}
