#!/usr/bin/perl -w

use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::SimpleAlign ;

open OUT1 ,">output.cobalt.prova.parsed";
my $in_seq  = Bio::SeqIO->new(-file =>"output.cobalt.prova" , '-format' => 'Fasta');
	
	while ( my $seq = $in_seq->next_seq() ) {
		my $id_seq =$seq->id();
		my $out_seq=$seq->seq();
		if ($id_seq =~ /^uaccno/){
			print  OUT1 ">",$id_seq,"\n",$out_seq,"\n";
		}
		#~ print ">",$id_seq,"\n",$out_seq,"\n";
		
		
	}
close OUT1;





