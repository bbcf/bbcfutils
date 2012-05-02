#!/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

our $title = "getRestEnzymeOccAndSeq.pl\tMarion Leleu\nVersion 0.01\n";
our $scriptname = "getRestEnzymeOccAndSeq.pl";

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

pod2usage(1) if (!($inputFile && $length && $outputFile && $fragmentFile));

if ($logFile) {
  open LOG,">$logFile" or die "Could not open $logFile\n";
  print LOG 
    "in $scriptname\ninput file=$inputFile\n"
      ."primary restriction site=$restSite\n"
        ."secondary restriction site=$restSiteSecond\nlength=$length\n"
          ."will create fragment file=$fragmentFile\n"
            ."and segment file=$outputFile\nlog file=$logFile\n"
              ."first restriction site=".$restSite
                ."second restriction site=".$restSiteSecond;
}

open IN, "<$inputFile" or die "Could not open $inputFile\n";
open OUT, ">$outputFile" or die "Could not open $outputFile\n";
open FRAG, ">$fragmentFile" or die "Could not open $fragmentFile\n";

print OUT join("\t",
               ("idOcc","idOccInChr","Chromosome","StartOcc","StartSeq (upstream)",
                "EndSeq (upstream)","Sequence (upstream)","StartSeq (downstream)",
                "EndSeq (downstream)","Sequence (downstream)",
                "indexOfrestSiteSecondOcc (in upstream seq)",
                "indexOfDpnOcc (in downstream seq)"))."\n";
print FRAG join("\t",
                ("idFrag","Chromosome","StartFrag","EndFrag","Sequence from Start",
                 "Start (start seq)","End (start seq)","Sequence from End",
                 "Start (end seq)","End (start seq)","indexOfrestSiteSecondOcc"))."\n";

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
my $lengthSecRestSite=($restSiteSecond eq "")?0:length($restSiteSecond);
my $fragStatus="FragIsValid";
#new header from genRep: >2746_NC_000081.5 gi|149301884|ref|NC_000081.5|NC_000081 Mus musculus chromosome 15, reference assembly (C57BL/6J)
while (<IN>) {
  chomp $_;
  if ($_ =~ />/) {
    $nHeader++;
    $header=$_;
    if ($nHeader==1) {
      $prevHeader=$header;
    }
    if ($nHeader>1) {
      $chr = ($prevHeader =~ />(\S+)\s/)?$1:"NA";
      if ($logFile) {
        print LOG "chromosome $chr\tsequence length=".length($seq)."\t";
      }
      $offset=0;
      my $curOcc = index($seq, $restSite, $offset);
      my $prevOcc = "";
      while ($curOcc != -1) {
        my $startSeq = substr($seq,$curOcc-$length,$length);
        my $endSeq =  substr($seq,$curOcc+$lengthRestSite,$length);
        print OUT join("\t",
                       ($i,$j,$chr,$curOcc,($curOcc-$length),$curOcc,$startSeq,
                        ($curOcc+$lengthRestSite),($curOcc+$lengthRestSite+$length),$endSeq))."\n";
        if ($prevOcc ne "") {
          my $allOccrestSiteSecondInFrag = &getAllOcc(substr($seq,$prevOcc,$curOcc-$prevOcc),$restSiteSecond);	
          my $status=&getStatus($curOcc-$prevOcc,$length,$restSiteSecond,$allOccrestSiteSecondInFrag);
          if ($allOccrestSiteSecondInFrag =~ /;/) #has at least one occ of 2nd rest site
            {
              my @allOccSecRS=split(";",$allOccrestSiteSecondInFrag);
              my $firstOcc=$allOccSecRS[0];
              my $lastOcc=$allOccSecRS[-1];
              print FRAG join("\t",
                              ($k,$chr,($prevOcc+$lengthRestSite),$curOcc,
                               substr($seq,$prevOcc+$firstOcc-$length,$length),
                               ($prevOcc+$firstOcc-$length),($prevOcc+$firstOcc),
                               substr($seq,$prevOcc+$lastOcc+$lengthSecRestSite,$length),
                               ($prevOcc+$lastOcc+$lengthSecRestSite),
                               ($prevOcc+$lastOcc+$lengthSecRestSite+$length),
                               $allOccrestSiteSecondInFrag,
                               ($curOcc-($prevOcc+$lengthRestSite)),$status))."\n";
              $k++;
            }
        } 
        $prevOcc = $curOcc;
				
        $offset = $curOcc + 1;
        $curOcc = index($seq, $restSite, $offset);
        $i++; $j++;
      }
      if ($logFile) {
        print LOG "$j occurrences of $restSite found in current chromsome ($i in total)\n";
      }
      $prevHeader=$header;
      $seq = "";
      $j=1; 
			
    }
  } else {
    $seq .= $_;
  }
}
$chr = ($prevHeader =~ />(\S+)\s/)?$1:"NA";
if ($logFile) {
  print LOG "chromosome $chr\tsequence length=".length($seq)."\t";
}
$offset=0;
my $curOcc = index($seq, $restSite, $offset);
my $prevOcc="";
while ($curOcc != -1) {
  my $startSeq = substr($seq,$curOcc-$length,$length);
  my $endSeq =  substr($seq,$curOcc,$length);
  print OUT join("\t",
                 ($i,$j,$chr,$curOcc,($curOcc-$length),$curOcc,$startSeq,
                  ($curOcc+$lengthRestSite),($curOcc+$lengthRestSite+$length),$endSeq))."\n";
  if ($prevOcc ne "") {
    my $allOccrestSiteSecondInFrag = &getAllOcc(substr($seq,$prevOcc,$curOcc-$prevOcc),$restSiteSecond);
    my $status=&getStatus($curOcc-$prevOcc,$length,$restSiteSecond,$allOccrestSiteSecondInFrag);
    if ($allOccrestSiteSecondInFrag =~ /;/) #has at least one occ of 2nd rest site
      {
        my @allOccSecRS=split(";",$allOccrestSiteSecondInFrag);
        my $firstOcc=$allOccSecRS[0];
        my $lastOcc=$allOccSecRS[-1];
        print FRAG join("\t",
                    ($k,$chr,($prevOcc+$lengthRestSite),$curOcc,
                     substr($seq,$prevOcc+$firstOcc-$length,$length),
                     ($prevOcc+$firstOcc-$length),($prevOcc+$firstOcc),
                     substr($seq,$prevOcc+$lastOcc+$lengthSecRestSite,$length),
                     ($prevOcc+$lastOcc+$lengthSecRestSite),
                     ($prevOcc+$lastOcc+$lengthSecRestSite+$length),
                     $allOccrestSiteSecondInFrag,($curOcc-($prevOcc+$lengthRestSite)),
                     $status))."\n";
        $k++;
      }
  }
  $prevOcc = $curOcc;

  $offset = $curOcc + 1;
  $curOcc = index($seq, $restSite, $offset);
  $i++; $j++;
}
if ($logFile) {
  print LOG "$j occurrences of $restSite found in current chromosome ($i in total)\n";
}


close IN;
close OUT;
close FRAG;
if ($logFile) {
  print LOG "result file=$outputFile\nDone\n";
  close LOG;
}

sub getAllOcc()
  {
    my $seqToSearch=shift;
    my $subseqToFind=shift;
    #print "search $seqToSearch\n";
    my $offset=0;
    my $res="";
    if ($subseqToFind eq "") {
      return "--";
    }

    my $n=0;
    my $curOcc = index($seqToSearch,$subseqToFind,$offset);
    while ($curOcc != -1) {
      $res .= $curOcc.";";
      $offset = $curOcc + 1;
      $curOcc = index($seqToSearch,$subseqToFind,$offset);
      $n++;
    }
    if ($res eq "") {
      $res="--";
    }
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
    foreach my $curOcc (@allOcc) {
      if ($curOcc eq "--") {
        $midFragStatus="midFragNotValid";last;
      } if ($curOcc>$segLength && $curOcc<$fragLength-$segLength) {
        $midFragStatus="midFragIsValid"; last;
      }
    }
    if ($midFragStatus =~ /midFragIsValid/ && $startSeqStatus =~ /StartSegOK/ && $endSeqStatus =~ /EndSegOK/) {
      $fragStatus="FragIsValid";
    } else {
      if ($startSeqStatus =~ /StartSegOK/) {
        if ($endSeqStatus =~ /EndSegOK/) {
          $fragStatus="FragIsNotValid";
        } else {
          $fragStatus="FragIsValid_badEndFrag";
        }
      } else {
        if ($endSeqStatus =~ /EndSegOK/) {
          $fragStatus="FragIsValid_badStartFrag";
        } else {
          $fragStatus="FragIsNotValid";
        }
      }
    }
    if ($restSiteSecond eq "") {
      $startSeqStatus="StartSeqOK";$endSeqStatus="EndSeqOK";$midFragStatus="midFragIsValid";$fragStatus="FragIsValid";
    }
    return "$startSeqStatus\t$endSeqStatus\t$midFragStatus\t$fragStatus";
  }
