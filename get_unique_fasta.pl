#!/usr/bin/perl -w

########## Script per obtindre sequencies fasta uniques. Fatsx_collpaser canvia el nom de la sequencia. Amb aquest script########
##          mantenim el nom de la sequencia
#################################################################################################################################

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
my @names;
my %params;
 
open OUT, ">temp";	 ## ULL enviem a un arxiu temporal!
##### Obrim el segon arxiu, que sera el que conte moltes seqs fasta ####
	my $fasta = Bio::SeqIO->new ( "-file"   => "_UTR.fasta",  "-format" => "Fasta" );
	while( my $inseq = $fasta->next_seq ) {
		my $id_seq = $inseq->id();
		my $sec = $inseq->seq();
		chomp($id_seq);
		unless(exists($params{$sec})) {
 			$params{$sec} += 1;
 			print OUT ">".$id_seq,"\n".$sec."\n";
 		} 				
				
 	}

}
close OUT;
print "DONE";



