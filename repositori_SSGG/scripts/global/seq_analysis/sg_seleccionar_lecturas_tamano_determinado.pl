#!/usr/bin/perl -w

######################

## Date: 24/9/12
## By: Sheila Zuniga

######################



use Bio::SeqIO;



sub usage
{
	print " USO: filtra lecturas de un fichero fasta atendiendo a su longitud \n";
	print " Input: \n";
	print " 	\$1: Fichero fasta \n";
	print " 	\$2: Longitud minima de lectura \n";
	
	print "	Output: Fichero fasta con las lecturas filtradas por el tamano definido \n";
	exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}




my $in = Bio::SeqIO->new (-file =>  $ARGV[0], -format => "Fasta");
my $seq_out = Bio::SeqIO->new(-format => "Fasta");

while(my $seq = $in->next_seq)
{
        if ($seq->length > $ARGV[1])
	{
		$seq_out->write_seq($seq);
	}
}

