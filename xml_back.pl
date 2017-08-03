#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;
use Bio::Search::Hit::GenericHit;

open OUT, ">cluster.hits";
my @files=<*>;
my @files2 = grep (/.tab/, @files);

#~ print @files2;

my $name1;
my $name2;
my $desc1;
my $desc2;

foreach my $blast_report(@files2){
# Get the report
	my $searchio = new Bio::SearchIO (-format => 'blastxml', -file=>$blast_report);
#~ print $searchio,"\n"; 

		while(my $result = $searchio->next_result){
			my $id=$result->query_name(); 
	#~ print $id,"\n"; ### contigs names
			my @hits=$result->hits;
	#~ print $hits[0],"\n";  #####imprimeix el best hit
	
				foreach(@hits){
		#~ my $obj = Bio::Search::Hit::GenericHit->new ();
		#~ print $obj,"\n";
		 #~ my $hit_names=$hit->name();
					$name1 = $hits[0]->name();
					$name2 = $hits[1]->name();
					$desc1= $hits[0]->description();
					$desc2= $hits[1]->description();
					print OUT $id."\t".$name1."\t".$desc1,"\n",$id."\t".$name2."\t".$desc2."\n";
				}	
		#~ print OUT $id."\t".$name."\t".$desc,"\n"; 
    #~ while(my $hits = $result->next_hit) {
		#~ print $hits,"\n"; ###### imprimeix tots els hits d'un contig
		
	   #~ my $name= $hits[0]->name();
	   #~ my $desc = $hits[0]->description();
	   #~ print $query."\t".$name."\t".$desc,"\n";
		#~ while (my $hsp=$hits->next_hsp){
		
#~ print  ">".$hit_name."\t".$desc."\t".$evalue."\n";
#~ print OUT ">$a $hit_name\n","$hit_string\n";
		}
	
	}
#~ }
close OUT;



