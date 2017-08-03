#!/usr/bin/perl  
use strict;

sub usage {
    print "\nPARA QUE SIRVE: Convierte las lecturas fastq a fasta\n";
    print "Input: 1) Fichero fastq  \n";
    print "Output: Fichero fasta   \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}


open(FASTA,"<",@ARGV[0]);
my $contador=0;
while  (my $fasta_linea=<FASTA>)
{
        chomp($fasta_linea);
	$contador = $contador+1;
	if ( $contador == 1 )
	{
		$fasta_linea=~s/@/>/;
		print  $fasta_linea, "\n";
	}
	elsif ( $contador == 2)
	{
		print  $fasta_linea, "\n";
	}
	elsif ($contador == 4) 
	{
		$contador=0;
	}

}
close (FASTA);
