#!/usr/bin/perl -w
# Calcula la media y la desviaci�n est�ndar de un intervalo conocido del archivo coverage.txt.
# El archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde a una posici�n
# y el valor es el correspondiente coverage de la posici�n
# El resultado aparece en pantalla.
# Se imprime un archivo con el nombre coordenadas_coverage.txt que contiene tres columnas. La primera nos da el valor
# de la coordenada en la secuencia total. La segunda, la posici�n en el intervalo seleccionado. En la tercera, 
# el valor del coverage para esa posici�n

# Si no introducimos el ning�n fichero, se imprime el modo de uso
sub usage {
	print "\n\nCOMO SE USA: sg_desviacion_estandar_coverage_de_exones.pl fichero\n\n";
	print "\nCalcula la media y la desviaci�n est�ndar de un intervalo conocido del archivo coverage.txt.";
    	print "\nEl archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde";
	print "\na una posici�n y el valor es el correspondiente coverage de la posici�n.\n";
 	print "El programa te pide la posici�n de inicio y la posici�n final del intervalo.\n\n";
	print "Si quieres calcular lo mismo de varios exones, puedes utilizar sg_desviacion_estandar_coverage_con_fichero_exones.pl\n\n";
	print "Se imprime el archivo coordenadas_coverage.txt, que contiene tres columnas. La primera nos da el valor de la coordenada\n";
	print "en la secuencia total. La segunda, la posici�n en el intervalo seleccionado. En la tercera, el valor\n";
	print "del coverage para esa posici�n\n\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

#Definimos las variables

# Pedimos los datos de inicio y final del intervalo del c�lculo

print "\nPrimera posicion: ";
my $inicio= <STDIN>;
print "\nUltima posicion: ";
my $final = <STDIN>;

# Abrimos el archivo
open(COVERAGE,"<",$ARGV[0]);

# Coge el intervalo e introduce cada posici�n en una posici�n del array @coverage
while  (my $linea=<COVERAGE>)
{
        chomp($linea);
        push (@coverage,$linea);
}

close(COVERAGE);

# Abrimos el archivo donde imprimiremos las coordenadas con el coverage
open(OUTPUT,">","coordenadas_coverage.txt") or die "No se puede abrir el archivo.\n";

# Imprimimos la cabecera del archivo de salida
print OUTPUT "coordenada_seq","\t","coordenada_intervalo","\t","coverage","\n";

# Definimos la variable $suma que nos dar� el resultado de la suma de las posiciones del array
my $suma = 0;
my $coordenada_intervalo = 1;

for ($i = $inicio - 1; $i <= $final - 1; $i++)
{
	$suma = $suma + $coverage[$i];
	my $coordenada_seq = $i + 1;
	print OUTPUT "$coordenada_seq","\t","$coordenada_intervalo","\t","$coverage[$i]","\n";
	$coordenada_intervalo++;
	$coordenada_seq++;
}

# Longitud del intevalo
$intervalo = $final - $inicio + 1;

# Calculamos la media
$media = $suma/$intervalo;

# Definimos la variable $cuadrado que nos dar� el resultado del cuadrado de la resta de la desviacion estandar
my $cuadrado = 0;

# Calculamos la desviacion estandar
for ($i = $inicio - 1; $i <= $final - 1; $i++)
{
	$parentesis = ($coverage[$i] - $media)**2;
	$cuadrado = $cuadrado + $parentesis;
}

# Calculamos la desviaci�n estandar
$raiz =(1/$intervalo)*$cuadrado;
$desviacion = sqrt ($raiz);

# Imprimimos los resultados
print "\nResultados\n";
print "Intervalo: ".$intervalo,"\n";
print "Media: ".$media,"\n";
print "Desviacion estandar: ".$desviacion,"\n\n";
print "El archivo coordenadas_coverage.txt contiene las coordenadas de tu secuencia con sus respectivos valores de coverage.\n\n";

# Cerramos el archivo de escritura
close (OUTPUT);
# Salimos del programa
exit;
