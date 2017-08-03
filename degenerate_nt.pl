#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

# Lookup table of degenerate IUPAC nucleotide codes.
my %deg2nuc = (
    "R" => ["A", "G"],
    "Y" => ["C", "T"],
    "S" => ["G", "C"],
    "W" => ["A", "T"],
    "K" => ["G", "T"],
    "M" => ["A", "C"],
    "B" => ["C", "G", "T"],
    "D" => ["A", "G", "T"],
    "H" => ["A", "C", "T"],
    "V" => ["A", "C", "G"],
    "N" => ["A", "C", "G", "T"]
);

# Recursive function that replaces degenerate nucleotides with all combinations.
sub generate
{
    if ($_[0] =~ /(.*)([RYSWKBDHVN])(.*)/) {
        #~ my $head = $1;
        #~ my $tail = $3;
        my @seqs;
        foreach my $nuc (@{$deg2nuc{$2}}) {
            #~ push @seqs, generate($head.$nuc.$tail);
	    push @seqs, $nuc;
        }
        return @seqs;
    }
    #~ else {
        #~ return $_[0];  ####### Si no hi ha cap nt degenerat retorna blanc
    #~ }
}

#Fasta file to be asked
print "Enter your fasta file name:\n";
 my $file = <STDIN>;
my $in  = Bio::SeqIO->new(-file => $file , '-format' => 'Fasta')  or die "Failed to open input file: $!";
   while ( my $seq = $in->next_seq() ) {
	     my $seq = $in->next_seq();
	             my $dna = $seq->seq();
		                       chomp($dna);

#Demo:  Print all sequences generated from ANCRG.
#~ print join("\n", generate("AACCG")), "\n";

# Print results:
#~ print join ("\n",generate ($dna));
print generate ($dna),"\n";
}
