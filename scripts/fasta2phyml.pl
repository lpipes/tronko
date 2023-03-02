#! /usr/bin/env perl

use strict;

foreach my $fname ( @ARGV ) {

  print "converting $fname...\n";
 my @spl = split(/\./,$fname);
  my $basename = `basename $fname .fasta`;
  chomp $basename;
  #my $outname = $basename.".phymlAln";
	my $outname = $spl[0] . ".phymlAln";
  open IN, "$fname" or die;
  open OUT, ">$outname" or die;

  sub pad {
    my ($char, $size) = @_;
    for(my $i=0; $i < $size; $i++) {
      print OUT "$char";
    }
  }

  my $maxSeqLen;
  my $maxNameLen;
  my $currSeqName = "";
  my $currSeq = "";
  my %seqHash;
  my @seqNames;

  foreach my $line ( <IN> ) {
    chomp $line;
    if($line =~ /^>/) {

      # call sub completeSeq (!)
      if($currSeqName eq "") {
	$currSeq eq "" or die;
      }
      else {
	$seqHash{$currSeqName} = $currSeq;
	push @seqNames, $currSeqName;
	if($maxSeqLen < length $currSeq) { $maxSeqLen = length $currSeq; }
	$currSeq = "";
	$currSeqName = "";
      }

      $line =~ s/^>//;
      if($maxNameLen < length $line) { $maxNameLen = length $line; }
      $currSeqName = $line;
    }
    else {
      $currSeq = $currSeq . $line;
    }
  }

  # call sub completeSeq (!)
  if($currSeqName eq "") {
    $currSeq eq "" or die;
  }
  else {
    $seqHash{$currSeqName} = $currSeq;
    push @seqNames, $currSeqName;
    if($maxSeqLen < length $currSeq) { $maxSeqLen = length $currSeq; }
    $currSeq = "";
    $currSeqName = "";
  }

  printf OUT "  %d %d\n", (scalar @seqNames), $maxSeqLen;
  foreach my $seqName ( @seqNames ) {
    print OUT $seqName;
    pad " ", (4 + $maxNameLen - (length $seqName));
    my $seqLine = $seqHash{$seqName};
    print OUT $seqLine;
    pad "-", ($maxSeqLen - (length $seqLine));
    print OUT "\n";
  }

  close IN;
  close OUT;

}

# completeSeq: for some reason the perl sub thing wasn't working like i wanted. 
#  sub completeSeq {
#    if($currSeqName eq "") {
#      $currSeq eq "" or die;
#    }
#    else {
#      print "pushed $currSeqName\n";
#      $seqHash{$currSeqName} = $currSeq;
#      push @seqNames, $currSeqName;
#      if($maxSeqLen < length $currSeq) { $maxSeqLen = length $currSeq; }
#      $currSeq = "";
#      $currSeqName = "";
#    }
#  };
