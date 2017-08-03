#!/usr/bin/perl  

sub usage {
    print "\nPARA QUE SIRVE: Convierte las lecturas truncadas de Fluxsimulator de longitud menor a la establecida anadiendo \"A\" hasta completar la longitud\n";
    print "EL ORDEN ES MUY IMPORTANTE \n";
    print "Input: 1) Fichero fastq; 2) longitud  \n";
    print "Output: Fichero fastq con las lecturas  \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(FASTA,"<",@ARGV[0]);
$longitud = @ARGV[1];
$contador=0;
while  (my $fasta_linea=<FASTA>)
{
        chomp($fasta_linea);
	$contador = $contador+1;
	if ( ($contador == 1) || ($contador == 3) )
	{
		print  $fasta_linea, "\n";
	}
	else 
	{
		$longlinea=length($fasta_linea);
		$nt_faltan= $longitud-$longlinea;
	
		if ($longlinea ne $longitud)
		{
			print $fasta_linea;
			if ($contador == 2)
			{
				for ($i=0; $i < $nt_faltan ; $i++)
                       		{
                               		print "A";
                       		}
			}
			elsif ($contador == 4)
			{
				$ultimo_QV=substr ($fasta_linea,$longlinea-1,1);
				for ($i=0; $i < $nt_faltan ; $i++)
				{
					print $ultimo_QV;
				}
				$contador=0;
			}	
			print "\n";	
		}
		
		else
		{
			print $fasta_linea, "\n";
			if ($contador == 4)
			{
				$contador=0;
			}
		}
	}

}





close (FASTA);
