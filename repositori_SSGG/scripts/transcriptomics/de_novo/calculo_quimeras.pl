#!/usr/bin/env perl

# -----------------------                                                                                                              
# -- Include libraries --
# -----------------------

use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;

# --------------------------------
# -- Input parameters & options --
# --------------------------------

my $lineA = $ARGV[0];
my $bam = $ARGV[1];

if (@ARGV != 2){
    die " 
    	Provide:
        	  	Input contig fasta file: fasta file to be renamed according to a specific size
        	  	Sample Bam file: from which contig name willbe retrieved
        	  	
\n".localtime()."\n=================================================================\n\n";
}

### Reads by bam file ###
my $length = system("samtools idxstats $bam > reads_$bam.count");


#### Parse contigs contigs file ### 
open IN, "<reads_$bam.count";
my @in =<IN>;
my %hash;

foreach my $line (@in){
    my @line_array = split ("\t" , $line);
    $hash{$line_array[0]} = $line_array[2]; #### Number of reads for each contig
    
}
#print Dumper(%hash);


my $input  = Bio::SeqIO->new(
			     -format => 'fasta',
			     -file   => $lineA,
			     );

while ( my $seq = $input->next_seq() ) {

	#open OUT ,">temp";
    my $sequence = $seq->seq;
    my $id = $seq->id;
    if ($hash{$id} != 0) {
	my $length = $hash{$id}; #### The number of reads
	#print OUT ">".$id."\n".$sequence,"\n";
	open OUT ,">$id.fasta";
	#system("formatdb -i temp -t $id.database -n $id.database -p F");
	#system("blastall -p blastn -d $id.database -i temp -e 10 -m 8 -o $id.output");
	print OUT ">".$id.";size=".$length.";\n".$sequence."\n";
    }
}
system ("cat Contig*.fasta > all_contigs.fasta");
system ("rm Contig*.fasta");


=head1 NAME
calculo_quimeras.pl

=head1 SYNOPSIS

calculo_quimeras.pl contigs.fasta sample.bam
