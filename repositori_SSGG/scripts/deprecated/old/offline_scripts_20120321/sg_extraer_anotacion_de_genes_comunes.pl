#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script extrae la informaci�n de un fichero de anotaci�n ya creado utilizando s�lo el nombre del gen. Est� pensado para extraer la informaci�n de la anotaci�n final de los genes comunes entre dos muestras.\n";
        print "\nCOMO SE USA: sg_extraer_anotacion_de_genes_comunes.pl <fichero_genes> <fichero_anotacion>\n";
        print "ejemplo: sg_extraer_anotacion_de_genes_comunes.pl genes_comunes.txt SNVs_out_SNP_effect.txt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero es la lista de genes. El segundo fichero es el de la anotaci�n.\n\n";
	print "OUTPUT: Fichero con tantas l�neas como l�neas donde aparezca el gen que buscamos. Hay que darle un fichero de salida.\n";
        exit(1);
}	

# Si s�lo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos los archivos de entrada
open(GENES,"<",$ARGV[0]);
open(ANOTACION,"<",$ARGV[1]);


# Cargamos en memoria el fichero de la anotaci�n del ensembl
for(my $j=1; my $anotacion=<ANOTACION>;$j++)
{
	chomp($anotacion);
	push(@array_anotacion,$anotacion);
}

# Cerramos el fichero del Ensembl
close(ANOTACION);

while($gen=<GENES>)
{
	chomp($gen);
	@split_gen = split(/\t/, $gen);
	for($n = 0; $n<=$#array_anotacion;$n++)
	{
		my @line_anotacion = split(/\t/,$array_anotacion[$n]);
		if($split_gen[0] eq $line_anotacion[10])
		{
			print $array_anotacion[$n],"\n";
		}
	}
}

exit;
