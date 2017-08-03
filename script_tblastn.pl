#!/usr/bin/perl -w
use strict;

my $a = "Waiting for a query name...\n";
	print $a;
my $line = <STDIN>;
chomp($line);
system("blastall -p tblastn -d $line.class.txt.database  -i ref.fasta -e 0.001 -m 7 -o ../../tblastn_results/$line ");
#print $line;
