#!/usr/bin/perl -w
# Calcula la media y la desviación estándar de los intervalos de secuencia que indiquemos en el segundo archivo
# conocido del archivo coverage.txt.
# El archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde a una posición
# y el valor es el correspondiente coverage de la posición
# El segundo archivo debe ser el archivo donde le indicamos los intervalos de secuencias y su nombre
# El resultado se imprime en el archivo exones_coverage.txt. La primera columna es el nombre de exón. La segunda,
# es la longitud. La tercera, la media y la cuarta, la desviación estándar. Para el cálculo, se han eliminado aquellos
# valores que estén fuera del intervalo media +- desviación estándar.

# Si no introducimos el ningún fichero, se imprime el modo de uso

sub usage {
	print "\n\nCOMO SE USA: sg_desviacion_estandar_coverage_con_fichero_exones.pl <fichero_coverage> <fichero_exones>\n";
	print "ejemplo: sg_desviacion_estandar_coverage_con_fichero_exones.pl coverage.txt exones.txt\n\n";
	print "Calcula la media y la desviación estándar del coverage de los exones utilizando el archivo coverage.txt.";
	print "Nos da los valores de media y desviación estándar eliminando los valores que se encuentren fuera del intervalo: media +- desviación estándar.";
    	print "\nEl primer archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponda";
	print "\na una posición y su correspondiente valor de coverage.\n";
	print "\nEl segundo archivo de entrada debe ser un fichero con tres columnas. La primera debe ser el número o nombre de exón.\n";
	print "El valor longitud indica la longitud de la secuencia que el programa resta a cada extremo de los intervalos.\n";
        print "La segunda, la coordenada de inicio del cada exón. La tercera, la coordenada del final.\n";
  	print "El programa imprime cinco columnas. La primera es el nombre o número de exon. La segunda, la longitud\n";
	print "del exón. La tercera, la media del coverage y la cuarta, la desviación estándar del coverage. La quinta, el porcentaje de datos excluídos del cálculo\n\n";
	print "Si sólamente quieres conocer estos datos de un intervalo, puedes utilizar sg_desviacion_estandar_coverage_de_exones.pl\n\n";
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

# Definimos la variable $suma que nos dará el resultado de la suma de las posiciones del array
my $suma;

# Definimos la variable $cuadrado que nos dará el resultado del cuadrado de la resta de la desviacion estandar
my $cuadrado;

print OUTPUT "exon","\t","longitud","\t","media coverage","\t","SD coverage","\t","% datos excluidos","\n"; 

# Calculamos la media y la desviación estándar para cada exón
for($j = 0; $j <= $#posiciones; $j++)
{
	$suma = 0;
	my $suma2 = 0;
	@ex = split (/\t/,$posiciones[$j]); # Cada columna de la fila de @posiciones la introducimos en una fila de @ex

        
# Calculamos la media y la desviación con todos los datos

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

	# Calculamos la desviación estandar
	$raiz =(1/$intervalo)*$cuadrado;
	$desviacion = sqrt ($raiz);
	$maximo = $media + $desviacion;
	$minimo = $media - $desviacion;

# Confirmamos los valores y volvemos a calcular la media y la desviación
	my $intervalo2 = 0; # El intervalo2 lo marcarán las veces que se ejecute el bucle if
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

	# Calculamos la desviación estándar
	
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
