############ Script per obtindre sequencies fasta (AA) de tots els hits obtinguts d'un resultat blastX ############
###################################################################################################################
 
#!/usr/bin/perl -w
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
print "Introdueix toxina:";

my $line = <STDIN>;
chomp($line);
my $searchio = new Bio::SearchIO (-format => 'blastxml', -file   => $line.".xml");
open OUTFILE , ">$line.fasta.txt" or die "can't open outfile\n"; ### Obrim el fitxer que servirà per tenir dues sequencies fasta: la de referencia i la problema

while( my $result = $searchio->next_result ) {
	#my $a=0;
	while(my $hit = $result->next_hit){
 		my $hit_name = $hit->name;
			#print $hit_name,"\n";
 			my $query = $result->query_name();
			#print $query,"\n";
		    my $desc=$hit->description();
				while (my $hsp = $hit->next_hsp ){
				    my $percentid = $hsp->percent_identity();
					my $score = $hit->bits();
					my $start = $hsp->start('hit');
					my $end = $hsp->end('hit');		
					my $hit_string = $hsp->hit_string;
					my $homo_string = $hsp->homology_string ;
					my $query_string = $hsp->query_string;
					my $frame = $hsp->frame('hit-frame');
					if ($hit_string =~ s/-//g ){   ####### Substitueixo els guions per tal que no em doni problemes a COBALT
						print OUTFILE ">".$hit_name."_".$frame."\n".$hit_string,"\n";
 					}else{
			    		print OUTFILE ">".$hit_name."_".$frame."\n".$hit_string,"\n";
   					}
						  			  			
				}		
		}
	}
}
print "DONE";
close OUTFILE;
