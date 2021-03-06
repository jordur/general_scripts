#!/usr/bin/perl

sub usage
{
	$usage=<<END;

- PARA QUE SIRVE: 
	Convierte los intervalos solapantes de un archivo BED de captura en �nicos.
	
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

for (my $j=0; $j<=$#fichero;$j++)
{
   @fila = split (/\t/,$fichero[$j]);
   @fila_anterior = split (/\t/,$fichero[$j - 1]);
   $resta = $fila_anterior[2] - $fila[1];
   $longitud_fila=$fila[2]-$fila[1];
   $valor=(abs($resta)/$longitud_fila)*100;
   if ($#fichero == 0)
   {
	
	print $fichero[$j],"\n";
   }
   else
   {
	if($j == 0)
	{
		print  $fila[0],"\t",$fila[1],"\t";
	}
	elsif($j == $#fichero)
	{
		if(($resta < -1) || (((abs($resta)/$longitud_fila)*100) < 40))
                {
                        print $fila_anterior[2],"\n",
			      $fila[0],"\t",$fila[1],"\t",$fila[2],"\n";
                }
		else
		{
			print  $fila[2],"\n",
		}
	}
	else
	{
		if(($resta < -1) || (((abs($resta)/$longitud_fila)*100) < 40))
		{
			
			print $fila_anterior[2],"\n",
		        $fila[0],"\t",$fila[1],"\t";
		}
	}
   }
}

exit;
