#!/usr/bin/perl -w

sub usage
{
	print "\n"; 
	print " PARA QUE SIRVE: Recorta las lecturas y sus valores de calidad asociados a la longitud especificada \n";
	print " COMO SE USA: \n";
	print " Input: 1) Fichero .csfasta o _QV.qual, uno de los dos, en formato FASTA 2) longitud que queremos que tengan las lecturas \n";
   	print " Output: el output es un fichero en formato FASTA. Hay que capturar la salida, es en pantalla\n\n";
    	exit(1);
}

if(scalar(@ARGV) == 0)
{
        usage();
}

open (CSFASTA, $ARGV[0]);
$longitud = $ARGV[1];

sub procesar_fastas
{
	$linea = shift;
	$primer_caracter = substr ($linea,0,1);
	if ($primer_caracter =~ m/(A|G|C|T)/)
	{
		$linea_recortada = substr ($linea,0,$longitud+1);
		print $linea_recortada,"\n";
	}
	else
	{
		@split_linea = split (" ",$linea);
		for ($i=0;$i<$longitud;$i++)
		{
			print $split_linea[$i]," ";
		}
		print "\n";
	}

}


sub main
{
	
	while ($lineas = <CSFASTA>)
	{
		chomp ($lineas);
		$primer_char=substr ($lineas,0,1);
		if ($primer_char ne "#")
		{
			if ($primer_char eq ">")
			{
				print $lineas,"\n";
			}
			else
			{
				&procesar_fastas($lineas);
			}

		}	
	}
}


&main();
close (CSFASTA);
