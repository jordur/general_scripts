#!/usr/bin/perl -w

use List::MoreUtils qw(uniq);

sub usage {
	print "\nEste script crea un archivo tabulado con la relación de genes únicos y el número de posiciones encontradas en cada gen.\n";
        print "\nCOMO SE USA: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl <fichero> \n";
        print "ejemplo: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl snp_no_descritos_chr2.txt\n\n";
	print "INPUT: fichero tabulado en el que la primera columna es el cromosoma, la segunda la coordenada de la SNV y la columna 11 es en la que se encuentra el nombre del gen. Está pensado para realizar el análisis de herencia recesiva o dominante con los ficheros de entrega de resecuenciación dirigida.\n\n";
	print "OUTPUT: fichero tabulado en el que la primera columna es el código del gen y la segunda, el número de variaciones encontradas por gen. Hay que darle un fichero de salida.\n";
        exit(1);
	}	
# Si sólo ejecutamos el script, se imprime las instrucciones de uso

if(scalar(@ARGV) == 0){
    usage();
}

# Abrimos el archivo 
open(NONSYN,"<",$ARGV[0]);

# Cargamos el fichero en memoria porque lo vamos a utilizar dos veces
for(my $f = 1; my $nonsyn=<NONSYN>; $f++)
{
	chomp($nonsyn);
	push(@no_sinonimas,$nonsyn);
}

# Cerramos el archivo de entrada
close(NONSYN);

# Generamos el array de genes que encontramos en las variaciones
for(my $j = 0; $j<=$#no_sinonimas; $j++)
{
	my @ex = split(/\t/,$no_sinonimas[$j]);
#	$posiciones = $ex[0]."-".$ex[1];
	if($ex[0] ne "Coordinates")
	{
		# Generamos un array con todos los genes
		push(@genes,$ex[10]);
	}
}
		# Esta variable nos va a dar el número de genes
		my $genes_unicos = uniq(@genes);
		# En este array tenemos el listado de genes no repetidos
		my @genes_unicos = uniq(@genes);

# Trabajamos con el array de genes no repetidos para contar el número de variaciones que encontramos en cada uno

for(my $i=0; $i<$genes_unicos; $i++)
{
	my @coordenadas;
	my $gen = $genes_unicos[$i];
	for(my $k = 0; $k<=$#no_sinonimas; $k++)
	{
        	my @ex2 = split(/\t/,$no_sinonimas[$k]);
#		my $chr_posiciones_gen = $ex2[0],"-",$ex2[1];
#		print $chr_posiciones_gen,"\n";
#		print "gen analizado en la posición $i...",$gen,".....","\n";
		if($ex2[10] eq $gen)
		{
#			print $ex2[10],"-------",$gen,"-----",$ex2[1],"\n";
			push(@coordenadas,$ex2[1]);
#			push(my @coordenadas,$chr_posiciones_gen);
#			print "variaciones en el array----",@variaciones,"\n";
		}
#		@coordenadas_unicas = uniq (@coordenadas);
	}
	my $variaciones_unicas = uniq (@coordenadas);
#	print "variaciones unicas---$variaciones_unicas\n";
	print "$gen\t$variaciones_unicas\n";

}

exit;
