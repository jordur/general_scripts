#!/usr/bin/perl

sub usage
{
	$usage=<<END;

- PARA QUE SIRVE: 
	Convierte los intervalos solapantes de un archivo BED de captura en únicos.
	
- COMO SE USA: sg_unificar_intervalos_de_captura_solapantes.pl sorted_chr1

- INPUT: fichero tabulado con los intervalos de captura. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera la final. Debe estar ordenado de menor a mayor.

ejemplo: chr1	343432	343445
	 chr1	343445	343645

- OUTPUT: fichero tabulado en el que la primera columna es el cromosoma, la segunda es el inicio del intervalo y, la tercera, es el final del intervalo\nejemplo: chr1	343432	343645

END
	print $usage;
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FILE,"<",@ARGV[0]);

my @fichero;
my @fila;
my @fila_anterior;

for (my $i=1; my $lines=<FILE>; $i++)
{
  chomp($lines);
  push (@fichero,$lines);
}
close (FILE);
my $valor=0;
for (my $j=0; $j<=$#fichero;$j++)
{
        @fila = split (/\t/,$fichero[$j]);
	@fila_anterior = split (/\t/,$fichero[$j - 1]);
	$resta = $valor - $fila[1];
#	print "\n -----> anterior $valor posterior $fila[1] resta $resta \n";
	#Si es la primera linea imprimo los valores de cromosoma y coordenada de inicio
	if($j eq 0)
	{
		print $fila[0],"\t",$fila[1],"\t";
		$valor = $fila[2];
	}

	#Final del fichero
	elsif($j == $#fichero)
	{
		#Si ultima linea no solapa con la anterior imprime la anterior y la ultima linea completa
		if($resta < -1)
                {
			print $valor,"\n",      
			$fila[0],"\t",$fila[1],"\t",$fila[2],"\n";
                }
		#Si la ultima linea solapa con la anterior imprimo el valor de la segunda columna de la ultima linea
		else
		{
			print  $valor,"\n",
		}
	}
	else
	{
		#El programa no entra aqui a no ser que las lineas no solapen, en tal caso imprime la segunda coordenada
		if($resta < -1)
		{
			print $valor,"\n",
		        $fila[0],"\t",$fila[1],"\t";
			$valor=$fila[2];
		}
		elsif (($resta > 0) && ($valor > $fila[2]))
		{
			$valor = $fila_anterior[2];
		}
		else
		{
			$valor = $fila[2];
		}
	}
}

exit;
