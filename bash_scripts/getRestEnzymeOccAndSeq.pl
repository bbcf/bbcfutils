#!/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

our $title = "getRestEnzymeOccAndSeq.pl\tMarion Leleu\nVersion 0.01\n";
our $scriptname = "getRestEnzymeOccAndSeq.pl";

#
#
#
#

my $inputFile;
my $restSite;
my $restSiteSecond="";
my $length;
my $outputFile;
my $fragmentFile;
my $logFile;

GetOptions(
        'i=s' => \$inputFile,
	'm=s' => \$restSite,
	's:s' => \$restSiteSecond,
	'l=i' => \$length,
        'o=s' => \$outputFile,
	'f=s' => \$fragmentFile,
	'x:s' => \$logFile,
) or pod2usage(2);

pod2usage(1) if (!$inputFile);
pod2usage(1) if (!$length);
pod2usage(1) if (!$outputFile);
pod2usage(1) if (!$fragmentFile);

if($logFile){open LOG,">$logFile" or die "Could not open $logFile\n";
print LOG "in $scriptname\n";
print LOG "input file=$inputFile\n";
print LOG "primary restriction site=$restSite\n";
print LOG "secondary restriction site=$restSiteSecond\n";
print LOG "length=$length\n";
print LOG "will create fragment file=$fragmentFile\n";
print LOG "and segment file=$outputFile\n";
print LOG "log file=$logFile\n";

print LOG "first restriction site=".$restSite;
print LOG "second restriction site=".$restSiteSecond;
}

open IN, "<$inputFile" or die "Could not open $inputFile\n";
open OUT, ">$outputFile" or die "Could not open $outputFile\n";
open FRAG, ">$fragmentFile" or die "Could not open $fragmentFile\n";

print OUT "idOcc\tidOccInChr\tChromosome\tStartOcc\tStartSeq (upstream)\tEndSeq (upstream)\tSequence (upstream)\tStartSeq (downstream)\tEndSeq (downstream)\tSequence (downstream)\tindexOfrestSiteSecondOcc (in upstream seq)\tindexOfDpnOcc (in downstream seq)\n";
print FRAG "idFrag\tChromosome\tStartFrag\tEndFrag\tSequence from Start\tStart (start seq)\tEnd (start seq)\tSequence from End\tStart (end seq)\tEnd (start seq)\tindexOfrestSiteSecondOcc\n";

my $header="";
my $prevHeader="";
my $nHeader=0;
my $seq="";
my $offset=0;
my $start;
my $chr;
my $end;
my $i=1;
my $j=1;
my $k=1;
my $lengthRestSite=length($restSite);
my $fragStatus="FragIsValid";
#>NC_000076.5 chr|NC_000076|NC_000076.5 Chromosome 10; [Mus musculus] gi:149288869 chromosome 10 complete sequence
#new header from genRep: >2746_NC_000081.5 gi|149301884|ref|NC_000081.5|NC_000081 Mus musculus chromosome 15, reference assembly (C57BL/6J)
while(<IN>)
{
	chomp $_;
	if($_ =~ />/){
		$nHeader++;
		$header=$_;
		if($nHeader==1){$prevHeader=$header;}
		if($nHeader>1)
		{
			$chr = ($prevHeader =~ /(.*)hromosome (.*),(.*)/)?$2:"chromosome not found";
			if($prevHeader =~ /mitochondrion/){$chr="chrM"}
			if($chr eq "chrMt"){$chr="chrM"}elsif($chr eq "Mt"){$chr="M"}
			if($logFile){print LOG "chromosome ".$chr."\t";}
			if($logFile){print LOG "sequence length=".length($seq)."\t";}
			$offset=0;
			my $curOcc = index($seq, $restSite, $offset);
			my $prevOcc = "";
			while ($curOcc != -1) {
               	 		my $startSeq = substr($seq,$curOcc-$length,$length);
		                my $endSeq =  substr($seq,$curOcc+$lengthRestSite,$length);
				print OUT $i."\t".$j."\t".$chr."\t".$curOcc."\t".($curOcc-$length)."\t".$curOcc."\t".$startSeq."\t".($curOcc+$lengthRestSite)."\t".($curOcc+$lengthRestSite+$length)."\t".$endSeq."\n";
				if($prevOcc ne "")
				{
					my $allOccrestSiteSecondInFrag = &getAllOcc(substr($seq,$prevOcc,$curOcc-$prevOcc),$restSiteSecond);
					my $status=&getStatus($curOcc-$prevOcc,$length,$restSiteSecond,$allOccrestSiteSecondInFrag);
					print FRAG $k."\t".$chr."\t".($prevOcc+$lengthRestSite)."\t".$curOcc."\t".substr($seq,$prevOcc+$lengthRestSite,$length)."\t".($prevOcc+$lengthRestSite)."\t".($prevOcc+$lengthRestSite+$length)."\t".substr($seq,$curOcc-$length,$length)."\t".($curOcc-$length)."\t".$curOcc."\t".$allOccrestSiteSecondInFrag."\t".($curOcc-($prevOcc+$lengthRestSite))."\t".$status."\n";
					$k++;
				} 
				$prevOcc = $curOcc;
				
				$offset = $curOcc + 1;
				$curOcc = index($seq, $restSite, $offset);
				$i++; $j++;
  			}
			if($logFile){print LOG $j." occurrences of ".$restSite." found in current chromsome ($i in total)\n";}
			$prevHeader=$header;
			$seq = "";
			$j=1; 
			
		}
	}
	else{$seq .= $_;}
}
$chr = ($prevHeader =~ /(.*)hromosome (.*),(.*)/)?$2:"chromosome not found";
if($prevHeader =~ /mitochondrion/){$chr="chrM"}
if($logFile){print LOG "chromosome".$chr."\t";
print LOG "sequence length=".length($seq)."\t";
}
$offset=0;
my $curOcc = index($seq, $restSite, $offset);
my $prevOcc="";
while ($curOcc != -1) {
	my $startSeq = substr($seq,$curOcc-$length,$length);
	my $endSeq =  substr($seq,$curOcc,$length);
	print OUT $i."\t".$j."\t".$chr."\t".$curOcc."\t".($curOcc-$length)."\t".$curOcc."\t".$startSeq."\t".($curOcc+$lengthRestSite)."\t".($curOcc+$lengthRestSite+$length)."\t".$endSeq."\n";
        if($prevOcc ne "")
        {
	 	my $allOccrestSiteSecondInFrag = &getAllOcc(substr($seq,$prevOcc,$curOcc-$prevOcc),$restSiteSecond);
		my $status=&getStatus($curOcc-$prevOcc,$length,$restSiteSecond,$allOccrestSiteSecondInFrag);
	        print FRAG $k."\t".$chr."\t".($prevOcc+$lengthRestSite)."\t".$curOcc."\t".substr($seq,$prevOcc+$lengthRestSite,$length)."\t".($prevOcc+$lengthRestSite)."\t".($prevOcc+$lengthRestSite+$length)."\t".substr($seq,$curOcc-$length,$length)."\t".($curOcc-$length)."\t".$curOcc."\t".$allOccrestSiteSecondInFrag."\t".($curOcc-($prevOcc+$lengthRestSite))."\t".$status."\n";		
		$k++;
        }
        $prevOcc = $curOcc;

	$offset = $curOcc + 1;
        $curOcc = index($seq, $restSite, $offset);
        $i++; $j++;
}
if($logFile){print LOG $j." occurrences of ".$restSite." found in current chromosome ($i in total)\n";}


close IN;
close OUT;
close FRAG;

if($logFile){
print LOG "result file=$outputFile\n";
print LOG "Done\n";
close LOG;
}

sub getAllOcc()
{
	my $seqToSearch=shift;
	my $subseqToFind=shift;
	#print "search $seqToSearch\n";
	my $offset=0;
	my $res="";
	if($subseqToFind eq ""){return "--";}

	my $n=0;
	my $curOcc = index($seqToSearch,$subseqToFind,$offset);
	while($curOcc != -1)
	{
		$res .= $curOcc.";";
		$offset = $curOcc + 1;
		$curOcc = index($seqToSearch,$subseqToFind,$offset);
		$n++;
	}
	#if($n>0){print $n." occurrences of restSiteSecond found\n"};
	if($res eq ""){$res="--"}
	return $res;

}

sub getStatus()
{
	my $fragLength=shift;
	my $segLength=shift;
	my $restSiteSecond=shift;
	my $allOccSecRestSite=shift;
	
	my @allOcc=split(";",$allOccSecRestSite);
	my $startSeqStatus = ($allOccSecRestSite eq "--" || $allOcc[0]>$segLength)?"StartSegOK":"StartSegNotValid";
        my $endSeqStatus = ($allOccSecRestSite eq "--" || $allOcc[$#allOcc]<$fragLength-$segLength)?"EndSegOK":"EndSegNotValid";
        my $fragStatus="FragIsValid";
	my $midFragStatus="midFragNotValid";
	foreach my $curOcc (@allOcc){if($curOcc eq "--"){$midFragStatus="midFragNotValid";last;} if($curOcc>$segLength && $curOcc<$fragLength-$segLength){$midFragStatus="midFragIsValid"; last;}}
#	if($midFragStatus =~ /NotValid/ || ($startSeqStatus =~ /NotValid/ && $endSeqStatus =~ /NotValid/)){$fragStatus="FragIsNotValid"}
	if($midFragStatus =~ /midFragIsValid/ && $startSeqStatus =~ /StartSegOK/ && $endSeqStatus =~ /EndSegOK/){$fragStatus="FragIsValid";}
	else
	{
		if($startSeqStatus =~ /StartSegOK/)
		{
			if($endSeqStatus =~ /EndSegOK/){$fragStatus="FragIsNotValid"}
			else{$fragStatus="FragIsValid_badEndFrag";}
		}
		else
		{
			if($endSeqStatus =~ /EndSegOK/){$fragStatus="FragIsValid_badStartFrag";}
			else{$fragStatus="FragIsNotValid"}
		}
	}
	if($restSiteSecond eq ""){$startSeqStatus="StartSeqOK";$endSeqStatus="EndSeqOK";$midFragStatus="midFragIsValid";$fragStatus="FragIsValid"}
	return $startSeqStatus."\t".$endSeqStatus."\t".$midFragStatus."\t".$fragStatus;
}
