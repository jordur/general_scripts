#!/usr/bin/perl -w
##### Script per parsejar blastxml agafant el primer hsp coma  best hit ####
## NOTA: l'escript es va fer dins el projecte de C.simus(mexican) que es troba a /media/bec2-jcalvete/Disc_Dur/Feina_Jordi/glandulas_mejicanas/glandulas_mejicanas/Feina/transcriptomes/global_results/clusters/CS/fsa_sequences/blastresults i el resultat bo el donava si el ID coincidia amb el HSP ID, d'aqui l'expressió regualr que trobem més avall 

use strict;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;
use Bio::Search::Hit::GenericHit;
use Bio::Search::HSP::HSPI;
use Bio::Search::Result::GenericResult;
use Bio::Search::HSP::GenericHSP;

open OUT, ">>Table_feautures.txt";
my @files=<*>;
my @files2 = grep (/.xml/, @files);

#print @files2;

my $start;
my $end;
my $strand;
my $a;
my $strand_name;

chomp($start);
chomp($end);
chomp($strand);
chomp($strand_name);

foreach my $blast_report(@files2){
# Get the report
	my $searchio = new Bio::SearchIO (-format => 'blastxml', -file=>$blast_report);
	while(my $result = $searchio->next_result){
	    my $id = $result->query_name(); 
	    (my $id2 = $id) =~ s/(.+)#.*/$1/;
	    chomp($id2);
#	    $end = $result->end('hit');
	    $a = 0;
#	    print $id2,"\n";
	    while(my $hit=$result->next_hit){
		$start = $hit->start('query');
		$end = $hit->end('query');
		if ($a < 1){		
		    while (my $hsp =$hit-> next_hsp){
			my $name1 = $hit->name();
 			(my $name2 = $name1 ) =~ s/\w+\|(\w+\d+)\|.+/$1/;  ### Expressió regualr per fer coincidir hspID amb sequenceID
			(my $name3 = $name1) =~ s/\w+\|(\w+\d+)\/.+/$1/; ### Idem anterior però en un altre format
#			print $id2,"\n";
#			print OUT $name1,"\n";
#			print $name2,"\n";
#			print $name3,"\n";
			if(($id2 =~ $name2)|| ($id2 =~ $name3)){
			    $a = $a + 1;
#			    print $name2,"\n";
			    $strand=$hsp->strand('query');
			    if ($strand == -1){
				$strand_name = "Minus";
			    }else{
				$strand_name = "Plus";
			    }
			print  OUT $id.":"."\tfrom"."\t".$start."\t"."to\t".$end."\t".$strand_name."\n";
#			print $id.":"."\tfrom"."\t".$start."\t"."to\t".$end."\n";
			}
		    }
		 }	
	    }
	}
}
close OUT;
