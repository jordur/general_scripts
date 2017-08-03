#!/usr/bin/perl -w

# Este script analiza las secuencias referencias y le da un valor de diez a las secuencias en may�sculas y un valor de cero a las bases en min�sculas. As� podremos calcular los intervalos de zonas repetitivas con el script sg_coverage.pl.
# Su principal utilidad es detectar las zonas repetitivas de las secuencias referencia, representadas por las letras min�sculas.
# El script se debe utilizar con las secuencias si validar.
# El resultado es un fichero de una columna en el que cada fila es una coordenada de la secuencia referencia, llamado zonas_repetitivas.txt
# Si no introducimos el ning�n fichero, se imprime el modo de uso.

sub usage {
        print "\nCOMO SE USA: sg_seleccionar_minusculas_y_mayusculas_de_secuencia_referencia.pl <secuencia_referencia> \n\n";
        print "ejemplo: sg_seleccionar_minusculas_y_mayusculas_de_secuencia_referencia.pl chr1.fasta\n\n";
        print "Analiza el fichero fasta y le da un valor de diez a las bases en may�sculas y un valor de cero a las bases en min�sculas.";
        print "Su principal utilidad es detectar las zonas repetitivas de las secuencias referencia, representadas por las letras min�sculas.";
        print "\nEl archivo de entrada debe ser un fichero fasta de la secuencia que quieres analizar. Hay que";
        print "\ntener en cuenta que no se deben utilizar las secuencias validadas.\n\n";
	print "El resultado es un fichero de una columna en el que cada fila es una coordenada de la secuencia referencia llamado zonas_repetitivas.txt\n\n";
	exit(1);
}

# Si no introducimos el fichero de entrada, se imprime el modo de uso
if(scalar(@ARGV) == 0){
    usage();
}

#Definimos las variables
my $secuencia_linea;
my $posicion_total = 1;

# Abrimos el archivo de entrada
open(SECUENCIA,"<",$ARGV[0]);

# Abrimos el archivo de salida
open (OUTPUT, ">", 'zonas_repetitivas.txt') or die "No puedo abrir el fichero de salida\n";

# Lee la secuencia l�nea a l�nea y no coge las l�neas que empiecen por ">"
while  (my $secuencia_linea=<SECUENCIA>)
{	
	my $primer_caracter = substr ($secuencia_linea,0,1);
	if ($primer_caracter ne ">")
        {
		# Eliminamos todos los espacios
		$secuencia_linea =~ s/\s//g;
		# Elimina el �ltimo espacio
		chomp($secuencia_linea);
		# Para cada l�nea de la secuencia, sustrae base a base. Si es min�scula, le asigna el valor cero,
		# si es may�sula, le asigna el valor diez
		for($i = 0 ; $i < length $secuencia_linea; ++$i)
		{
			$base = substr($secuencia_linea, $i, 1);
			if ($base =~ m/[a-z]/)
			{
				$id = 0;
			}
			else
			{
				$id = 10;
			}
			# Imprime cada valor en una l�nea diferente
			print OUTPUT $id,"\n";
			$posicion_total++;
		}
	}
}

# Cerramos el archivo de entrada
close(SECUENCIA);

# Cerramos el archivo de salida
close(OUTPUT);

#Salimos del programa
exit;
