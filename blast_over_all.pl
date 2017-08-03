#!/usr/bin/perl -w

use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;
#####Script per blastejar diferentes seq entre elles##########
############Obtenim les nt fasta d'una secuencia donada############
my @list=qw/16124242 28972958 62131097 46451422 186659794 3426323 88861961 6492259 31742524 14915765 62463 167614288/;
# 
# 
 foreach my $prot (@list){
 open OUT, ">$prot.fa";
 my $out = Bio::SeqIO->new ( "-file"   => ">prot.fa",
 			    "-format" => "Fasta" );
 my $gb = Bio::DB::GenBank->new();
# 
 	my $gi = <STDIN>;
 	chomp($gi);
 # 	
 	print "Fetching gi '", $prot, "' ... ";
	my $seq = $gb->get_Seq_by_id($prot); # GI Number
 	print $seq,"\n";
 	my $seqio = $seq->seq();
 	print OUT ">",$prot,"\n",$seqio,"\n\n";
# 
 close OUT;
		}
###########Creem database de cada secuencia nt de la llista#############
 system("for file in `ls *.fa`; do formatdb -i $file -p F -n $file.database;done");


####### Blast de cada sequencia fasta contra una db determinada (en cada cas la db=la propia sequencia) #############3
 system ("for file in `ls *.fa`; do blastall -p blastn -d 16124242.fa.database -i $file -e 0.001 -m 9 >> all_over_161 ; done");
