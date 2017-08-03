#!/usr/bin/perl
use List::Util qw (sum);

#Inputs:
#coverage.txt (1 columna) ó fwd_rev_coverage.txt (2 columnas)
#Outputs:
#Nos interesa establecer las posiciones de inicio y final de las subcadenas mayores del triple de la media (p.e. 232), asi como el tamano de 
#cada una de estas subcadenas.

my $numero_bases;#Calcular el numero de bases de la secuencia de referencia.
my @fila;#Del coverage.txt.
my @suma;#De cada fila anterior.
my @posiciones;# De las delecciones.
my @inicio;#Array que contendra las coordenadas iniciales de cada deleccion.
my @final;#Array que contendra las coordenadas finales de cada deleccion.
my @longitud;#Array que contiene el tamaño de las deleccione.


print "\nCoverage minimo: ";
my $coverage_maximo= <STDIN>;
print "\nNumero minimo de bases entre areas: ";
my $numero_bases_entre_areas = <STDIN>;
print "\nCalculando high coverage ...\n";


#Suma al array SUMA.

#open (FILE1, "<", 'coverage.txt');
open (FILE1, "<", @ARGV);#Fichero de entrada como parametro. I/O en ruta local usuario.

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  push(@suma,@fila[0]+@fila[1]);
}

close (FILE1);#Cierro la entrada del fichero coverage.txt


#Coger las posiciones de las lineas del array SUMA que son 0 y las guarda en array POSICIONES.

for (my $j=0;$j<=$#suma;$j++)
{
  if ($suma[$j] >= $coverage_maximo)#Solo cojo aquellas lineas del fichero mayores de (media*3=232) 
  {
    push(@posiciones,$j+1);
  }
  $numero_bases=$j+1;
}                                                               

open (OUT, ">", 'out_high_coverage_areas.xls') or die "No puedo abrir out_high_coverage_areas.xls.\n";
print OUT "Inicio\tFinal\tTamano\n";

if ($#posiciones!=-1)
{

  #En el array de indices POSICIONES separo los indices de inicio (a INICIO) y final (a FINAL) de cadenas-de-columna con todo "0". 

  @inicio=(@posiciones[0]);#El primer valor del array INICIO tiene que ser el primero de POSICIONES y como nunca va a cumplir el if lo anado de forma manual

  for (my $i=0;$i<=$#posiciones+1;$i++)
  {	
        #Solo aquellos numeros en los que la diferencia entre ese numero y el siguiente sea mayor o igual a 2	
	if ((@posiciones[$i+1]-@posiciones[$i])>=$numero_bases_entre_areas) 
	{	
                push(@final,@posiciones[$i]);#Guardo las coordenadas del final de la delecion
		push(@inicio,@posiciones[$i+1]);#Guardo las coordenadas del inicio de la delecion
	}    
  }

  push(@final,@posiciones[-1]);#El ultimo valor del array FINAL tiene que ser el ultimo de POSICIONES, de otra forma nunca lo obtendria
	

  #Imprime el principio, el final y el tamano de la deleccion en un fichero llamado output.xls (y el tamano de nuevo en el array LONGITUD).

  my $tamano;#Contiene el tamano de la deleccion

  for (my $a=0;$a<=$#inicio;$a++)
  {	
	$tamano=(@final[$a]-@inicio[$a])+1;#Calcula el tamano de la delecion
	print OUT @inicio[$a]."\t".@final[$a]."\t".$tamano."\n"; 
	push (@longitud,$tamano);#meto los tamanos de las delecciones en el array
  }	

  #print  OUT "El numero total de bases que no cambian  ", $numero_bases-sum(@longitud)," pb","\n";

}

print "\nOperacion finalizada.\n\n";

close (OUT);#Cierro out_high_coverage_areas.xls
