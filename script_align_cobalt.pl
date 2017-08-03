#!/usr/bin/perl -w
######################################################
# Script per alinear totes les seqüències de toxines #
# fasta obtingudes a una proteïna de referència.     #
######################################################

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

###### Demanem el directori de treball com a argument ################
if (@ARGV < 1){
	print "Afegir com a argument el directori de treball\n"
}else{
	
#### Fixem el directori de treball com l'argument #########
chdir $ARGV[0];


#######Esborrem referencies anteriors si existeixen ############

if (-e "output.cobalt.prova"){
	system("rm output.cobalt.prova");
}elsif(-e "output.cobalt.prova.parsed"){
	system ("rm output.cobalt.prova.parsed");
}

####### introducir protein query and reference id per tal d'obtenir des de NCBI el CDS de la proteïna i guardar-lo a ref.fasta ########
my $prot_id;
my $sec_ref;
my $id;	    	    

print "Waiting for a query name...:\n";
my $line = <STDIN>;	
chomp($line);
my $gb = Bio::DB::GenBank->new();

print "Waiting for reference id...:\n";
my $gi = <STDIN>;
chomp($gi); 


##### Si existeix ref.fasta no faig res, d'altra forma la busco a NCBI ##############
if (-e "$line.ref.fasta"){
	print "Existeix $line.ref.fasta. Anem a BLAST\n";
}else{

	print "Fetching gi ", $gi, " ... ";
	my $seq = $gb->get_Seq_by_gi($gi); # GI Number
	$sec_ref = $seq->seq();
	$id = $seq->display_id();
	print "done."."\n";
			    
	my @features = $seq->get_SeqFeatures(); # just top level
	foreach my $feat ( @features ) {
		if( $feat->primary_tag eq 'CDS' ) {
			( $prot_id )= $feat->get_tag_values('protein_id');		    
			print "Fetching prot id ", $prot_id, " ... ";
			my $prot_obj = $gb->get_Seq_by_acc( $prot_id );
		#	my $prot_obj = $gb->get_Seq_by_acc( $gi );
			print "done."."\n";
			print "Writting fasta ...";
  	    	$sec_ref=$prot_obj->seq();
 	    	print "done."."\n";

			open REF,">$line.ref.fasta";
        	print REF ">".$id."\n";
        	print REF $sec_ref."\n";
        	print "done\n";
#			     	     chomp(my $ref=<REF>);
#			     	     print REF ">".$prot_id."\n";
#			     	     print REF $sec_ref."\n";		 
			#print OUTFILE ">".$prot_id."\n";
			#print OUTFILE "$sec_ref\n";
        	close REF;

		}
				 
	}
}

####### tblastn . Si existeix el fitxer xml no faig res #######
if (-e "$line.xml"){
	print "Existeix el fitxer BLAST\n";
}else{
	print "Performing BLAST process..";
	#system("blastall -p tblastn -v 12000 -d database/${line}.database -i ref.fasta -e 0.001 -m 7 -o ${line}.xml");
	system("blastall -p tblastn -d database/${line}.database -i ${line}.ref.fasta -e 0.001 -m 7 -o ${line}.xml");
	print "done\n";
}


#######parseo blastxml####### 
chomp($line);
my $blast_report = "$line.xml";
chomp($blast_report);
#print $blast_report;

### Per si vull introduir el report manualment ###
#print "Enter blast_report file name\n";
#my $blast_report = <STDIN>;

my $searchio = new Bio::SearchIO (-format => 'blastxml',
 							      -file   => $blast_report);
#my $nt_seqs = new Bio::SeqIO (-format => 'fasta',
#							  -file   => 'pla2_all_correct_frame.fasta');
my $query;
my $qstart;
my $qend;
my $hstart;
my $hend;

open ORF , ">$line.orf" or die "can't open outfile\n"; ### Obrim el fitxer q contindrà els contigs NOMES amb les regions que donen BLAST positiu. A aquest contigs alinearem els reads per quantificar.

while( my $result = $searchio->next_result ) {
	#my $a=0;
    while(my $hit = $result->next_hit){
 		my $hit_name = $hit->name;
			#print $hit_name,"\n";
 			my $query = $result->query_name();
			#print $query,"\n";
		    my $desc=$hit->description();
				while (my $hsp = $hit->next_hsp ){
					open OUTFILE , ">$line.fasta.txt" or die "can't open outfile\n"; ### Obrim el fitxer que servirà per tenir dues sequencies fasta: la de referencia i la problema
				    my $percentid = $hsp->percent_identity();
					my $score = $hit->bits();
					my $start = $hsp->start('hit');
					my $end = $hsp->end('hit');		
					my $hit_string = $hsp->hit_string;
					my $homo_string = $hsp->homology_string ;
					my $query_string = $hsp->query_string;
					my $frame = $hsp->frame('hit-frame');
					$qstart = $hsp->start('query');
					$qend = $hsp-> end ('query');	
					$hstart = $hsp -> start ('hit');
					$hend = $hsp -> end ('hit');
					#$a++;
					#print $hit_name."\n".$hit_string,"\n";
					if ($hit_string =~ s/-//g ){   ####### Substitueixo els guions per tal que no em doni problemes a COBALT
						print OUTFILE ">".$hit_name."_".$frame."\n".$hit_string,"\n";
 					}else{
			    		print OUTFILE ">".$hit_name."_".$frame."\n".$hit_string,"\n";
   					}
					######### Obtindre la com a primera linea la referencia per tal que el COBALT sigui només un alinemanent 1 a 1. Fem l'alineament amb COBALT ############
					my $out;
					my $id;
					my $in  = Bio::SeqIO->new(-file => "$line.ref.fasta" , -format => 'Fasta');
					while ( my $seq = $in->next_seq() ) {
						my $id= $seq->id();
						my $out=$seq->seq();   
						print OUTFILE ">".$id."\n".$out."\n";
						print "Aligning with COBALT...$hit_name\n";
						#system("/home/bec2-jcalvete/Feina_Jordi/programes/bin/cobalt -rps_evalue 0.1 -rpsdb ~/Feina_Jordi/programes/cobalt/cdd_clique_0.75 -p ~/Feina_Jordi/programes/cobalt/patterns -blast_evalue 0.1 -max_dist 0.95 -i $line.fasta.txt >> ${line}_temp");
			  		} 	
					
					##### Amb això obtinc les posicions de començament i final del contig obtinguda del BLAST per tal d'obtindre els fitxers fasta .orf que em ###
					#### permetran alinear els reads només a aquelles sequencies codificants. Aixo es especielament relevant en el cas de les BPP ### 
					
					#my $dif = $hend - $hstart;					
					my $input = Bio::SeqIO->new(-file => "$line.fasta" , -format => 'Fasta');
					while ( my $seq = $input->next_seq() ) {
						my $sequencia = $seq->seq();
						my $id = $seq->id();
						if ($id eq $hit_name){
							print ORF ">".$id."\n".$seq->subseq($hstart,$hend)."\n";
							#print ">".$id."\n".$hstart,"\t".$hend."\n";
						}
					
					}
						  			  			
				}		
	}
}
close OUTFILE;
close ORF;
######## Sequence alignment ##########
#print "Aligning with COBALT...";
#system("/home/bec2-jcalvete/Feina_Jordi/programes/bin/cobalt -rps_evalue 0.1 -rpsdb ~/Feina_Jordi/programes/cobalt/cdd_clique_0.75 -p ~/Feina_Jordi/programes/cobalt/patterns -blast_evalue 0.1 -max_dist 0.95 -i $line.fasta.txt > ${line}_cobalt");
#print "done.\n";
#########

# print "Aligning with COBALT... ";
# my $seq_in= Bio::SeqIO->new ( "-file"   => "outfile.txt",  "-format" => "Fasta" );
# my $a = 1;
#      while( my $inseq = $seq_in->next_seq ) {
# 	 my $id = $inseq->id();
# 	 my $sec = $inseq->seq();
# 	 open IN, "<ref.fasta";
# 	 my @ref = <IN>;
# 	 open OUT, ">temp";
# 	 print OUT @ref;
# 	 print OUT ">".$id."\n".$sec."\n";
# 	 close OUT;
# 	 close IN;
# 	 system("./cobalt -rpsdb ~/programes/cobalt/cdd_clique_0.75 -clusters F -p ~/programes/cobalt/patterns -i temp >> temp2");
# 	 system("rm temp");
# }
# print "done.\n";


####### Parsejem la sortida de COBALT per obtindre un fitxer amb una sola referencia i tots els alineaments baix. Per aixo obtenim la sequencia de referencia de nou i aquelles que NO son iguals a la id de referencia ##########  
#open OUT, ">".$line."_cobalt";
my $id_ref;
my $id_seq;
my $ref = Bio::SeqIO->new ( "-file"   => "$line.ref.fasta",  "-format" => "Fasta" );
while( my $refseq = $ref->next_seq ) {
	$id_ref = $refseq->id();
	my $sec_ref = $refseq->seq();
	#print OUT ">".$id_ref."\n".$sec_ref."\n";
}
#print $id_ref,"\n";
my $in= Bio::SeqIO->new ( "-file"   => "${line}_temp",  "-format" => "Fasta" );
while( my $inseq = $in->next_seq ) {
	$id_seq = $inseq->id();
   	if($id_seq ne $id_ref){
 		my $sec = $inseq->seq();
 		#print ">".$id_seq."\n".$sec."\n";
 		#print OUT ">".$id_seq."\n".$sec."\n";
     }
}
#system("rm ${line}_temp");
#close OUT;	
}
print "DONE\n";
