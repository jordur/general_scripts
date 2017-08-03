#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: Crear un fichero con valores de calidad a partir del nombre y longitud de un contig\n";
    print " Input: 1)Fichero con dos columnas, el nombre del contig y la longitud del mismo \n\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FASTA,"<",@ARGV[0]);
#open(MIRBASE,"<",@ARGV[1]);
#my $nombre_chr=@ARGV[2];
#open(SENSE,">","$nombre_chr.conocidos.sense.ma.tab") or die "No puedo abrir $nombre_chr.conocidos.sense.ma.tab.\n";
while (my $linea =<FASTA>)
{
	chomp ($linea);
	@split_linea = split ("\t",$linea);
	print ">".$split_linea[0],"\n";
	$contador=1;
	for (my $i=1; $i <= $split_linea[1]; $i++)
	{
		if ($i < $split_linea[1])
		{
			print "20 ";
		}
		elsif ($i= $split_linea[1])
		{
			print "20\n";
		} 
	}
	#$contador=1;
}
