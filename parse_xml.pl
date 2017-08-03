#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SearchIO;
use Bio::SeqFeature::SimilarityPair;
use Bio::Search::Hit::GenericHit;
use Cwd;
### Script per parsejar xml BLAST (outfmt 6) a TAB BLAST (outfmt 5). 
if (@ARGV < 1){
	print "Afegir com a argument el directori de treball"
	
}else{

### Obrim l'arxiu ouput al directori fixat com a argument 1 ###
chdir $ARGV[0];

my @files=<*>; ## Llegim tots els arxius 
my @files2 = grep (/^blastN_results_*/, @files); ## LLegim resultats blastx
#my @files2 = grep (/^blast_5_contigsX_*/, @files); ## LLegim resultats blastx

my $name;
my $desc;
my $desc_2;
my $acc;
my $evalue;
my $length;
my $bits;
my $query_start;
my $query_end;
my $subject_start;
my $subject_end;
my $percentid;
my $gaps;

foreach my $blast_report(@files2){

	$blast_report =~ /blastN_results_(\d+)/;
	#print $1,"\n";
	open OUT, ">>blast_output".$1;
######## Get the report
	my $searchio = new Bio::SearchIO (-format => 'blastxml', -file=>$blast_report);
#	print $searchio,"\n"; 
		while(my $result = $searchio->next_result){
			my $id=$result->query_name(); 
			$acc = $result->query_description();
			my $b =1;
#			print OUT $id,"\n"; ### contigs names
			while (my $hit=$result->next_hit){
			    if($b < 2){
					$name = $hit->name(); ### hit name
					$desc=$hit->description(); ### description name
					($desc_2 = $desc) =~ s/\s/_/g; ## Canvio els espais per underscores per no tenir problemes amb la delimitació de camps de Trinotate
					while( my $hsp = $hit->next_hsp ) {
						$evalue = $hsp ->evalue(); ### evalue
						$length = $hsp ->length('total'); #### alignment lenght
						$bits = $hsp -> bits(); ### bit score
						$query_start = $hsp->start('query'); ### query start
						$query_end= $hsp -> end ('query'); ### query end 
						$subject_start = $hsp->start('hit'); ### subject start
						$subject_end= $hsp -> end ('hit'); ### subject end
						$percentid = $hsp->percent_identity(); ### percentage identity
						$gaps = $hsp->gaps('total'); ### number of total gaps
						
					}
					
#				print $id."\t".$name."\t".$desc,"\n";
				print OUT $acc."\t".$desc_2."\t".$percentid."\t".$length."\t"."10"."\t".$gaps."\t".$query_start."\t".$query_end."\t".$subject_start."\t".$subject_end."\t".$evalue."\t".$bits."\n";
				
			    }
			    $b++;
				}	
		
		    }
		    
		}
close OUT;
print "Done"
}


