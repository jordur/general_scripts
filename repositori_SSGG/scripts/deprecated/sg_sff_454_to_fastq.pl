#!/usr/bin/perl

sub usage {
	my $usage =<<END;
- The fasta filename must end in either .fasta or .fna

- The quality filename must have the same basename as the fasta file and end with .qual. For example, if your sequence file is "foo.fna" then the quality file must be named "foo.qual".
- Link: http://seqanswers.com/forums/showthread.php?t=2775
END

        print $usage;
}

use warnings;
use strict;
use Data::Dumper;
use File::Basename;

if(scalar(@ARGV) == 0){
    usage();
}

my $inFasta = $ARGV[0];
my $baseName = basename($inFasta, qw/.fasta .fna/);
my $inQual = $baseName . ".qual";
my $outFastq = $baseName . ".fastq";

my %seqs;

$/ = ">";

open (FASTA, "<$inFasta");
my $junk = (<FASTA>);

while (my $frecord = <FASTA>) {
	chomp $frecord;
	my ($fdef, @seqLines) = split /\n/, $frecord;
	my $seq = join '', @seqLines;
	$seqs{$fdef} = $seq;
}

close FASTA;

open (QUAL, "<$inQual");
$junk = <QUAL>;
open (FASTQ, ">$outFastq");

while (my $qrecord = <QUAL>) {
	chomp $qrecord;
	my ($qdef, @qualLines) = split /\n/, $qrecord;
	my $qualString = join ' ', @qualLines;
	$qualString =~ s/\s+/ /g;
	my @quals = split / /, $qualString;
	print FASTQ "@","$qdef\n";
	print FASTQ "$seqs{$qdef}\n";
	print FASTQ "+\n";
	foreach my $qual (@quals) {
		print FASTQ chr($qual + 33);
	}
	print FASTQ "\n";
}

close QUAL;
close FASTQ;
