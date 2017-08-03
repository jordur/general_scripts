#!/usr/bin/perl -w                                                                            

use Bio::SeqIO;
use Bio::SearchIO;

print "Print a number\n";
my $counter = <STDIN>;
chomp($counter);
print "Print contig name\n";
my $c= <STDIN>;
chomp ($c);
#open OUT, ">".$c."_".$counter.".fasta";
print "Print start position\n";
my $a = <STDIN>;
chomp($a); 
print "Print end position\n";
my $b = <STDIN>;
chomp($b);
$in  = Bio::SeqIO->new(-file => $c , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) {
#    if ($seq->id ){
    print ">",$seq->id,"\n",$seq->subseq($a,$b),"\n";

#}
}
