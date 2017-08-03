#!/usr/bin/perl -w
### Aquest script chekeja que el comando split a la consola efectivament ha dividit els fastes un per linia. Puc fer split amb 
### split -l LINIES ARXIU_FASTA o amb sg_split... i donava resultats diferents. Amb això miro que cada arxiu resultant té un id i un seq.
### Arguments: DIRECTORI 

use Bio::SeqIO; 

open OUT, ">",$ARGV[0]."/errors.txt";
chdir $ARGV[0];

my @files=<*>; ## Llegim tots els arxius 
my @files2 = grep (/transfuse_25_cons.fa(\d+).fasta/, @files); ## LLegim fitxers fasta amb REGEX específica 

foreach my $file(@files2){
	#my $in = Bio::SeqIO->new (-file =>  @ARGV, -format => "fasta"); 
	my $in = Bio::SeqIO->new (-file =>  $file, -format => "fasta"); 
	while(my $seq = $in->next_seq){ 
		my $iden = $seq->primary_id;
		my $seq = $seq -> seq;
		if (defined $iden and defined $seq){
			next;
		}else{
			print OUT "Existeix un ERROR a l'arxiu ".$file."\n";
		} 
		#if ($iden =~ /Transcript_1$/){
			#print ">".$iden,"\n";
			#print $seq,"\n"; 
		#} 

	} 
}
print "Done";
#close OUT;