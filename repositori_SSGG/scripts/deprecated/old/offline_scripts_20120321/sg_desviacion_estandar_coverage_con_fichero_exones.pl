#!/usr/bin/perl -w
# Calcula la media y la desviaci�n est�ndar de los intervalos de secuencia que indiquemos en el segundo archivo
# conocido del archivo coverage.txt.
# El archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde a una posici�n
# y el valor es el correspondiente coverage de la posici�n
# El segundo archivo debe ser el archivo donde le indicamos los intervalos de secuencias y su nombre
# El resultado se imprime en el archivo exones_coverage.txt. La primera columna es el nombre de ex�n. La segunda,
# es la longitud. La tercera, la media y la cuarta, la desviaci�n est�ndar. Para el c�lculo, se han eliminado aquellos
# valores que est�n fuera del intervalo media +- desviaci�n est�ndar.

# Si no introducimos el ning�n fichero, se imprime el modo de uso

sub usage {
	print "\n\nCOMO SE USA: sg_desviacion_estandar_coverage_con_fichero_exones.pl <fichero_coverage> <fichero_exones>\n";
	print "ejemplo: sg_desviacion_estandar_coverage_con_fichero_exones.pl coverage.txt exones.txt\n\n";
	print "Calcula la media y la desviaci�n est�ndar del coverage de los exones utilizando el archivo coverage.txt.";
	print "Nos da los valores de media y desviaci�n est�ndar eliminando los valores que se encuentren fuera del intervalo: media +- desviaci�n est�ndar.";
    	print "\nEl primer archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponda";
	print "\na una posici�n y su correspondiente valor de coverage.\n";
	print "\nEl segundo archivo de entrada debe ser un fichero con tres columnas. La primera debe ser el n�mero o nombre de ex�n.\n";
	print "El valor longitud indica la longitud de la secuencia que el programa resta a cada extremo de los intervalos.\n";
        print "La segunda, la coordenada de inicio del cada ex�n. La tercera, la coordenada del final.\n";
  	print "El programa imprime cinco columnas. La primera es el nombre o n�mero de exon. La segunda, la longitud\n";
	print "del ex�n. La tercera, la media del coverage y la cuarta, la desviaci�n est�ndar del coverage. La quinta, el porcentaje de datos exclu�dos del c�lculo\n\n";
	print "Si s�lamente quieres conocer estos datos de un intervalo, puedes utilizar sg_desviacion_estandar_coverage_de_exones.pl\n\n";
	print "El resultado se imprime en el fichero exones_coverage.txt.\n\n";
	exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

# Abrimos el archivo coverage.txt
open(COVERAGE,"<",$ARGV[0]);

while (my $linea=<COVERAGE>)
{
        chomp($linea);
        push (@coverage,$linea);
}

# Abrimos el archivo con los exones y sus coordenadas de inicio y final
open(EXONES,"<",$ARGV[1]);

while (my $exon=<EXONES>)
{
        chomp($exon);
        push (@posiciones,$exon);
}

close(COVERAGE);

close(EXONES);

# Abrimos el archivo donde queremos imprimir los resultados
open(OUTPUT,">","exones_coverage.txt");

# Definimos la variable $suma que nos dar� el resultado de la suma de las posiciones del array
my $suma;

# Definimos la variable $cuadrado que nos dar� el resultado del cuadrado de la resta de la desviacion estandar
my $cuadrado;

print OUTPUT "exon","\t","longitud","\t","media coverage","\t","SD coverage","\t","% datos excluidos","\n"; 

# Calculamos la media y la desviaci�n est�ndar para cada ex�n
for($j = 0; $j <= $#posiciones; $j++)
{
	$suma = 0;
	my $suma2 = 0;
	@ex = split (/\t/,$posiciones[$j]); # Cada columna de la fila de @posiciones la introducimos en una fila de @ex

        
# Calculamos la media y la desviaci�n con todos los datos

	for ($i = $ex[1] - 1; $i <= $ex[2] - 1; $i++)
        {
                $suma = $suma + $coverage[$i];
        }
        # Longitud del intevalo
        $intervalo = $ex[2] - $ex[1] + 1;

	# Calculamos la media
	$media = $suma/$intervalo;

	$cuadrado = 0;

	# Calculamos la desviacion estandar
	for ($k = $ex[1] - 1; $k <= $ex[2] - 1; $k++)
	{
		$parentesis = ($coverage[$k] - $media)**2;
		$cuadrado = $cuadrado + $parentesis;
	}

	# Calculamos la desviaci�n estandar
	$raiz =(1/$intervalo)*$cuadrado;
	$desviacion = sqrt ($raiz);
	$maximo = $media + $desviacion;
	$minimo = $media - $desviacion;

# Confirmamos los valores y volvemos a calcular la media y la desviaci�n
	my $intervalo2 = 0; # El intervalo2 lo marcar�n las veces que se ejecute el bucle if
	for($h = $ex[1] - 1; $h <= $ex[2] - 1; $h++)	
	{
		if ( $coverage[$h] >= $minimo && $coverage[$h] <= $maximo)
		{
			$intervalo2 = $intervalo2 + 1;
			$suma2 = $suma2 + $coverage[$h];

        		# Calculamos la media
			$media2 = $suma2/$intervalo2;

		}
	}

	# Calculamos la desviaci�n est�ndar
	
	my $intervalo3 = 0;
	for($n = $ex[1] - 1; $n <= $ex[2] - 1; $n++)
        {
		$cuadrado2 = 0;
                if ( $coverage[$n] >= $minimo && $coverage[$n] <= $maximo)
                {
               		$intervalo3 = $intervalo3 + 1;
		        $parentesis2 = ($coverage[$n] - $media2)**2;
			$cuadrado2 = $cuadrado2 + $parentesis2;
			$raiz2 =(1/$intervalo3)*$cuadrado2;
			$desviacion_final = sqrt ($raiz2);

		# Calculamos el porcentaje de datos que excluimos
		$porcentaje = 100 - ($intervalo3/$intervalo)*100;


		}
	}
# $porcentaje = 100 - ($intervalo3/$intervalo)*100;
# Imprimimos los resultados
 print OUTPUT "$ex[0]","\t","$intervalo3","\t","$media2","\t","$desviacion_final","\t","$porcentaje","\n";
}

# Cerramos el fichero de salida
close(OUTPUT);

# Salimos del programa
exit;
