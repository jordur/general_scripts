#!/usr/bin/perl -w
# Calcula la media y la desviaci�n est�ndar de un intervalo conocido del archivo pileup.txt.
# El archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde a una posici�n
# y el valor es el correspondiente pileup de la posici�n
# El resultado aparece en pantalla.
# Se imprime un archivo con el nombre coordenadas_pileup.txt que contiene tres columnas. La primera nos da el valor
# de la coordenada en la secuencia total. La segunda, la posici�n en el intervalo seleccionado. En la tercera, 
# el valor del pileup para esa posici�n

# Si no introducimos el ning�n fichero, se imprime el modo de uso
sub usage {
	print "\n\nCOMO SE USA: sg_desviacion_estandar_pileup_de_exones.pl fichero\n\n";
	print "\nCalcula la media y la desviaci�n est�ndar de un intervalo conocido del archivo pileup.txt.";
    	print "\nEl archivo de entrada debe ser un fichero con una sola columna donde cada fila corresponde";
	print "\na una posici�n y el valor es el correspondiente pileup de la posici�n.\n";
 	print "El programa te pide la posici�n de inicio y la posici�n final del intervalo.\n\n";
	print "Si quieres calcular lo mismo de varios exones, puedes utilizar sg_desviacion_estandar_pileup_con_fichero_exones.pl\n\n";
	print "Se imprime el archivo coordenadas_pileup.txt, que contiene tres columnas. La primera nos da el valor de la coordenada\n";
	print "en la secuencia total. La segunda, la posici�n en el intervalo seleccionado. En la tercera, el valor\n";
	print "del pileup para esa posici�n\n\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

#Definimos las variables

# Abrimos el archivo
open(PILEUP,"<",$ARGV[0]);

# Coge el intervalo e introduce cada posici�n en una posici�n del array @pileup
while  (my $linea=<PILEUP>)
{
        chomp($linea);
        push (@pileup,$linea);
}

close(PILEUP);


# Definimos la variable $suma que nos dar� el resultado de la suma de las posiciones del array
my $suma = 0;

for ($i = 0; $i <= $#pileup; $i++)
{
	@fila = split (/\t/,$pileup[$i]);
	$suma = $suma + $fila[3];
}

$media = $suma/$i;
print "\nN�mero de bases: $i\nMedia: $media\n\n";


# Salimos del programa
exit;
