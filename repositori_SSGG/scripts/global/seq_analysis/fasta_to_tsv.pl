#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO; 
use Bio::SearchIO;
use File::Basename;
use Getopt::Long;
use Cwd;

my $file = $ARGV[0];
my $specie = $ARGV[1];
open OUT, ">".$specie.".tsv";


(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);


GetOptions (
                    "c=s"   => \$file,
                    "v=s"   => \$specie,
	    );

print "\n=================================================================
$basename: Template script for developing Perl scripts v1.1\n";

#if (not defined $file || not defined $specie) {
    die "
          Options:
          -c contig file
          -v Venomics_ID
\n".localtime()."\n=================================================================\n\n" unless defined $file && defined $specie;
#}

 
my $in  = Bio::SeqIO->new(-file => $file,
			  -format => 'Fasta');

print OUT "GLAND_ID\tSEQUENCE_ID\tSEQUENCE\tLENGTH\tSPECIE\n";
while ( my $seq = $in->next_seq() ) {
    my $id = $seq ->id;
    my $sequence = $seq -> seq;
    my $desc = $seq -> desc;
    my $length = $seq -> length;
    print OUT $specie."\t".$id."\t".$sequence."\t".$length."\t".$desc,"\n";
}
