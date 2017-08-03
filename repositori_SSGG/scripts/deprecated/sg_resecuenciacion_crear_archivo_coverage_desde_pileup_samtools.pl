#!/usr/bin/perl -w

# Script para crear los archivos de coverage a partir de archivos generados con el pileup de Samtools, en la que la cuarta segunda columna corresponde a la posición en la secuencia referencia y la cuarta al coverage. Crea dos ficheros de salida.
# 1. coverage.txt: compuesto por una columna en la que cada fila corresponde a una posición de la secuencia referencia.
# 2. posicion_coverage.txt: archivo tabulado de dos columnas. La primera contiene las posiciones de la secuencia referencia y la segunda, el coverage de cada posición.
# El archivo de entrada es el archivo de salida del pileup de Samtools y se introduce manualmente la longitud de la secuencia referencia.

sub usage {
	print "\nEste script está pensado para obtener los valores de coverage de cada posición en la secuencia referencia. Funciona correctamente al aplicarle la función pileup del paquete Samtools a los archivos sam generados en el mapeo con Bioscope.\n";
        print "\nCOMO SE USA: sg_resecuenciacio_crear_archivo_coverage_desde_pileup_samtools.pl <fichero_pileup_samtools> \n";
        print "ejemplo: sg_desviacion_estandar_coverage_con_fichero_exones.pl pileup_samtools.txt\n\n";
        print "Script para crear los archivos de coverage a partir de archivos generados con el pileup de Samtools, en la que la segunda columna corresponde a la posición en la secuencia referencia y la cuarta al coverage. Crea dos ficheros de salida:\n\n";
	print "1. coverage.txt: compuesto por una columna en la que cada fila corresponde a una posición de la secuencia referencia.\n\n";
	print "2. posicion_coverage.txt: archivo tabulado de dos columnas. La primera contiene las posiciones de la secuencia referencia y la segunda, el coverage de cada posición.\n\n";
	print "El archivo de entrada es el archivo de salida del pileup de Samtools sin la opción -c.\n Se introduce manualmente la longitud de la secuencia referencia.\n\n";
        exit(1);

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
}
if(scalar(@ARGV) == 0){
    usage();
}


# Pedimos la longitud de la secuencia referencia
print "\nLongitud de la secuencia referencia: ";
$longitud_referencia = <STDIN>;

# Abrimos el archivo coverage
open(COVERAGE,"<",$ARGV[0]);

while (my $coverage=<COVERAGE>)
{
        chomp($coverage);
        push (@posicion,$coverage);
}

close(COVERAGE);

# Abrimos el archivo de salida posicion_coverage.txt
open(OUTPUT,">","posicion_coverage.txt");

# Abrimos el archivo de salida coverage.txt
open(COV,">","coverage.txt");

my $k = 1;
for ($i= 0; $i <= $#posicion; $i++)
{	
	@ex = split(/\t/,$posicion[$i]);
	# Si la coordenada se encuentra en el archivo pileup, la imprime y pasa a la siguiente
	if( $k eq $ex[1])
	{
		print OUTPUT "$ex[1]","\t","$ex[3]","\n";
		print COV "$ex[3]","\n";
		$k++;
        }
	# Si no la encuentra, imprime los valores de las coordenadas y le asigna un valor cero de coverage
        else
        {
		for($j=$k; $j < $ex[1]; $j++)
		{
			print OUTPUT "$j","\t","0","\n";
			print COV "0","\n";
		}
		# Para pasar correctamente a la siguiente fila del archivo pileup, debemos imprimir la fila en la que estamos
		print OUTPUT "$ex[1]","\t","$ex[3]","\n";
		print COV "$ex[3]","\n";
		# Marcamos el valor de $k para que continúe en el valor de pileup siguiente al último impreso
		$k=$ex[1] + 1;	
	}
}

# Si las lecturas no mapean hasta el final de la referencia, tenemos que completar los archivos con valores de coverage iguales a cero, una vez ha acabado de extraer las filas del fichero pileup
for ($h = $k; $h <= $longitud_referencia; $h++)
	{
		print OUTPUT "$h","\t","0","\n";
		print COV "0","\n";
	}

# Cerramos los ficheros de salida
close(OUTPUT);
close(COV);

exit;

