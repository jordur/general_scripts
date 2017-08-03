#!/usr/bin/env perl

use strict;
use warnings;

my ( $file ) = @ARGV;

open INFILE, $file or die "can't open input:$!\n";
open OUTFILE, ">" . $file . "_output.txt" or die "can't open output:$!\n";

my %seen_pair;
while ( my $line = <INFILE> ){
  chomp $line;
  my ($col_1, $col_2) = split (/\t/, $line);
  if (  $seen_pair{ $col_1 }{$col_2} || $seen_pair{ $col_2 }{$col_1} ){
    print STDERR "Pair already seen: $col_1\t$col_2\n";
    next;
  }
  $seen_pair{ $col_1 }{$col_2} = 1;
  print OUTFILE "$col_1\t$col_2\n";
}

close INFILE;
close OUTFILE;


__END__

=head1 NAME - nrPair

=head1 AUTHOR

Javier Santoyo Lopez B<email> jsantoyo@cipf.es

=head1 DATE

September 16, 2010
