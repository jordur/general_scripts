#!/usr/bin/perl -w

use Bio::SeqIO; 
my $in = Bio::SeqIO->new (-file =>  @ARGV, -format => "Fasta"); 
while(my $seq = $in->next_seq) 
{ 
	print $seq->primary_id, "\t", $seq->length, "\n"; 
	#print $seq->seq, "\n"; 
} 
