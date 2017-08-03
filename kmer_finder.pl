#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

#Sort out input
die "Usage: $0 FASTAfile\n" if(!@ARGV);
my $filename = $ARGV[0];

my $kmer = 7;
my $lowerLimit = 0;

my %dataHash;


open(IN, "gzip -dcf $filename | ") or die "Could not open input file: $!\n";

while (my $line = <IN>) {
	chomp($line);
	next if($line eq "");
	next if($line =~ m/^>/);
	next if($line =~ m/^\+/);
	next if($line =~ m/^@/);
	next if($line !~ m/^[ACTGN]+$/i);
	
	my $len = length($line);
	#print $line,"\t",$len,"\n"; # imprimeix la sequencia nucleotídica i la longitud.
	for (my $i = 0; $i <= ($len-$kmer); $i++) {

		my $sub = substr($line, $i, $kmer);
		#print $sub,"\n"; ## va sumant una lletra més a la dreta
		#print STDERR " "x$i;
		#print STDERR $sub."\n";

		if(defined($dataHash{$sub})) {
			$dataHash{$sub} += 1;
		}
		else {
			$dataHash{$sub} = 1;
		}
	}
	#print STDERR "$line\n";
}
close IN;

#my $num = 0;
foreach my $j (sort {$dataHash{$a} <=> $dataHash{$b}} (keys(%dataHash))) {
	if($dataHash{$j} > $lowerLimit) {
#		print ">$num\n",$j,"\n";
#		$num = $num +1;
		print "$j\n";
#		print "$j\t$dataHash{$j}\n";
	}
}


