#!/usr/bin/perl -w
use File::ReadBackwards;

sub usage {
        print "\nCOMO SE USA: \n\n";
	print "INPUT: 1) fichero muestras tabulado, posicion y coverage 2) fichero con las posiciones de referencia\n";
	print "OUTPUT: columna con el coverage por posicion, si una base no se cubre el programa imprime un 0\n\n";
        exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}


#Cargo en memoria las posiciones de la referencia.
open(POSICIONES,$ARGV[$#ARGV]);
	
	while (my $posiciones=<POSICIONES>)
	{
		chomp($posiciones);
		push (@posiciones_referencia,$posiciones);
	}

close(POSICIONES);


$str = "comp";
#Abre el listado de ficheros correspondientes para generar la matriz, van separados por cromosoma y por muestra
for($j = 0; $j <= $#ARGV-1; $j++)
{
	
	open (SALIDA,">",$ARGV[$j].$str);
	open (FICHERO,"<",$ARGV[$j]) or die ("no encuentro el fichero - $!");
	#Veo la ultima linea del fichero en caso de que no apareciera nada en esa muestra en ese cromosoma de manera que el coverage correspondiente seria 0.
	my $ultima_linea_fichero  = File::ReadBackwards->new($ARGV[$j])->readline;
	chomp ($ultima_linea_fichero);
	my @split_ultima = split (/\s+/,$ultima_linea_fichero);
	if (!$ultima_linea_fichero)
	{
		for ($p=0; $p <= $#posiciones_referencia ; $p++)
		{
			print SALIDA "0\n";
		}	
	}
	else
	{
		my $val=0;
		while (my $lineas = <FICHERO>)
		{
			chomp ($lineas);
			@posiciones_muestra = split (/\s+/,$lineas);
			for ($mo=$val; $mo <= $#posiciones_referencia ; $mo++)
			{
				chomp ($posiciones_referencia[$mo]);
				#Solo hay dos posibilidades, o bien la posicion de la columna de referencia y la muestra son iguales o bien la posicion de la columna de la muestra es mayor que la columna de referencia, nunca puede ser mas pequenya porque las posiciones de las muestras son exclusivamente aquellas que estan en rango.
				if($posiciones_referencia[$mo] eq $posiciones_muestra[0])
				{
					#Si la posicion en la columna de referencia es igual a la ultima fila del fichero de entrada y el numero de posiciones que han de cubrirse es mayor entonces significa que hay bases que no se han cubierto y se imprime un coverage de 0.
					if (($mo < $#posiciones_referencia) && ($posiciones_referencia[$mo] eq $split_ultima[0]))
					{
						#Imprimo 1 porque la posicion en la columna de referencia es igual a la ultima fila
						print SALIDA $posiciones_muestra[1],"\n";
						#Para imprimir el resto de posiciones de la columna de referencia que no se cubren en el fichero de entrada resto al numero de posiciones en la columna de referencia el numero de lineas que ya he recorrido.
						$resto=$#posiciones_referencia-$mo;
						for ($quedan =0; $quedan <=$resto-1; $quedan++)
						{
							print SALIDA "0\n";
						}
						last;
					
					}
					else
					{
						print SALIDA $posiciones_muestra[1],"\n";
						$val=$mo+1;
						last;
					}
				}
				elsif ($posiciones_muestra[0] > $posiciones_referencia[$mo])
				{
					print SALIDA "0\n";
				}
			}
		}
	}
	close(FICHERO);
	close(SALIDA);
}










