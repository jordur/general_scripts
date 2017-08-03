#!/usr/bin/perl -w

########## Script per obtindre els ID d'un fitxer fasta, localitzar-los en un altre i obtindre la seva sequencia fasta ########
###############################################################################################################################

use strict;
use Bio::DB::GenBank;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::SimpleAlign ;

if (@ARGV < 1){
	print "Afegir com a argument el directori de treball\n"
}else{

chdir $ARGV[0];
my $id_ref;    
print "Waiting for a query name...:\n";
my $line = <STDIN>;	
chomp($line);
my @names;
my %params;
open OUT, ">".$line."_cobalt.fasta"; ####### Script on enviare les seqs fasta positives

##### Obrim el primer arxiu per obtindre els IDs i posar-los en un array #########
my $cobalt = Bio::SeqIO->new ( "-file"   => "${line}_cobalt",  "-format" => "Fasta" );
while( my $cobseq = $cobalt->next_seq ) {
	my $cobalt_seq = $cobseq->id();
	(my $cob_seq2 = $cobalt_seq) =~ s/_\d+$//g;
	chomp($cob_seq2);
	push @names, $cob_seq2; 
	#print $cob_seq2,"\n";
	%params = map { $_ => 1 } @names;
}
	
#### Obrim el segon arxiu, que sera el que conte moltes seqs fasta ####
	my $fasta = Bio::SeqIO->new ( "-file"   => "${line}.fasta",  "-format" => "Fasta" );
	while( my $inseq = $fasta->next_seq ) {
		my $id_seq = $inseq->id();
		#print $id_seq,"\n";
		if(exists($params{$id_seq})) {
 			my $sec = $inseq->seq();
 			#print OUT $id_seq,"\n";
 			#print OUT ">".$id_seq."\n".$sec."\n";
 			print  OUT ">".$id_seq."\n".$sec."\n";
    	 }
	}
	### Loop a través dels id de cobalt ####
#	foreach(@names){
#		print $_,"\n";
#	}	
}
print "DONE";
