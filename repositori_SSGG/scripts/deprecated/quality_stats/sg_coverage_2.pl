#!/usr/bin/perl -w

# Input:
# Archivo tabulado de una columna en que cada línea corresponde a cada posición de la secuencia referencia. El programa necesita dos valores. El primero, llamado "coverage máximo", y corresponde a el valor de coverage máximo que vamos a aceptar que esté incluído en una deleción. El segundo valor, llamado "Numero de bases entre deleciones", nos va a indicar la longitud mínima de una deleción para ser tenida en cuenta.
# Output:
# Un archivo tabulado, llamado out_coverage.txt en el que vienen detalladas las coordenadas de inicio, las coordenadas finales y la longitud de las deleciones

sub usage {
        print "\nCOMO SE USA: sg_coverage_2.pl <fichero_coverage> \n\n";
        print "ejemplo: sg_coverage_2.pl coverage.txt\n\n";
        print "Input:\nArchivo tabulado de una columna en que cada línea corresponde a cada posición de la secuencia referencia. El programa necesita dos valores. El primero, llamado 'coverage máximo', y corresponde al valor de coverage máximo que vamos a aceptar para que esté incluído en un intervalo de bajo coverage. El segundo valor, llamado 'Numero de bases entre deleciones', nos va a indicar la longitud mínima de los intervalos de bajo coverage para ser tenidos en cuenta.\n\n";
        print "Output:";
        print "\nUn archivo tabulado, llamado out_coverage.txt en el que vienen detalladas las coordenadas de inicio, las coordenadas finales y la longitud de los intervalos de bajo coverage";
        print "\n\nEjemplo:\n";
	print "Inicio  Final   Tamano\n";
	print "123     223     100\n\n";
        exit(1);
}

# Si no introducimos el fichero de entrada, se imprime el modo de uso
if(scalar(@ARGV) == 0){
    usage();
}


# Definimos el valor de coverage máximo para que acepte los valores de coverage
print "\nCoverage maximo [0 por defecto]: ";
$coverage_maximo = <STDIN>;
chop ($coverage_maximo);

if ($coverage_maximo le 0)
{
  $coverage_maximo = 0;
}

# Definimos el tamaño mínimo de las deleciones que queremos detectar
print "\nNumero de bases entre deleciones [2 por defecto]: ";
$num_bases_entre_delecciones = <STDIN>;
chop ($num_bases_entre_delecciones);

if ($num_bases_entre_delecciones le 2)
{
  $num_bases_entre_delecciones = 2;
}

# Calculamos la longitud de la secuencia referencia
$count = `wc -l < $ARGV[0]`;
die "wc failed: $?" if $?;
chomp($count);
$longitud_mas_uno = $count + 1;

print "Longitud de la secuencia referencia: ",$count,"\n";

# Abrimos el fichero de entrada y lo convertimos en un fichero temporal para que imprimir en el fichero temporal tmp_coverage.txt '0' cuando el valor de coverage sea menor que $coverage_maximo y '1' cuando sea mayor

# Abrimos el archivo temporal
open(TEMP, ">",'tmp_coverage.txt');

# Abrimos el archivo de entrada
open(INPUT, "<", $ARGV[0]);

print "Analizando los datos...";

# Convertimos los valores de coverage a nuestros valores relativos
while(my $linea=<INPUT>)
{ 
	chomp($linea);
	if($linea > $coverage_maximo)
	{
		$coverage_falso	= 1;
	}
	else
	{
		$coverage_falso = 0;
	}
	print TEMP $coverage_falso,"\n";
}

# Cerramos ambos ficheros
close (INPUT);
close(TEMP);

# Miramos cuál es el último valor de coverage 
my $ultima_base = system("tail -n 1 $ARGV[0]");
if($ultima_base > $coverage_maximo)
{
        $coverage_ultima = 1;
}
else
{
        $coverage_ultima = 0;
}

# Imprimimos el archivo con las deleciones. Para ello, analizamos el archivo temporal.
my $posicion = 1;

# Empezamos definiendo arbitrariamente el valor inicial de la variable $valor como '0'
my $valor = 0;
my $intervalo = 0;

# Abrimos el fichero de salida del programa e imprimimos la cabecera de las columnas
open(OUTPUT, ">",'out_coverage.txt');
print OUTPUT "Inicio","\t","Final","\t","Tamano","\n";

# Abrimos el archivo temporal y lo leemos línea a línea
open(FINAL,"tmp_coverage.txt");
while($lines=<FINAL>)
{
	chomp($lines);
	
	# Si nuestro valor de coverage es igual que el valor, incrementamos el intervalo
	if($valor eq $lines)
	{
		$intervalo++;
	}
	else
	{
		# Si la línea tiene el valor de '1' y el valor es '0', acaba la deleción, luego imprimimos los valores del intervalo
		if($valor == 0)
		{
			if($intervalo > $num_bases_entre_delecciones)
			{
				print OUTPUT $posicion - $intervalo,"\t",$posicion - 1,"\t",$intervalo,"\n";
			}
		$valor = 1;
		$intervalo=1;
		}
		else		# Si es el valor es '1' y la línea 0, empieza la deleción, luego tenemos que darle un valor de '0' a $valor para que empiece a aumentar el intervalo
		{
			$valor = 0;
			$intervalo=1;
		}
	}
	$posicion++;
	# Si los valores de nuestro coverage acaban en cero, debemos imprimirlos a parte, ya que se quedan fuera del bucle
	if(($longitud_mas_uno eq $posicion) && ($coverage_ultima eq 0) && ($intervalo > $num_bases_entre_delecciones))
	{
		print OUTPUT $posicion - $intervalo,"\t",$posicion - 1,"\t",$intervalo,"\n";
	}
}
close(FINAL);

print "\nOperación finalizada.\n\n";

system("rm tmp_coverage.txt");

# Salimos del programa
exit();
