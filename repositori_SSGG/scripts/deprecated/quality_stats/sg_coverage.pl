#!/usr/bin/perl
use List::Util qw (sum);
#use Smart::Comments '###';

#Inputs:
#coverage.txt (1 columna) ó fwd_rev_coverage.txt (2 columnas)
#Outputs:
#Las posiciones de inicio y final de las subcadenas con "0", asi como el tamano de cada una de estas subcadenas.

my $numero_bases;#Calcular el numero de bases de la secuencia de referencia.
my @fila;#Del coverage.txt.
my @suma;#De cada fila anterior.
my @posiciones;# De las delecciones.
my @inicio;#Array que contendra las coordenadas iniciales de cada deleccion.
my @final;#Array que contendra las coordenadas finales de cada deleccion.
my @longitud;#Array que contiene el tamaño de las deleccione.

print "\nCoverage maximo [0 por defecto]: ";
$coverage_maximo = <STDIN>;
chop ($coverage_maximo);

if ($coverage_maximo <= 0)
{
  $coverage_maximo = 0;
}

print "\nNumero de bases entre delecciones [2 por defecto]: ";
$num_bases_entre_delecciones = <STDIN>;
chop ($num_bases_entre_delecciones);

if ($num_bases_entre_delecciones <= 2)
{
  $num_bases_entre_delecciones = 2;
}
     
print "\nCargando el fichero  <".$ARGV[0]."> en memoria ...\n";
#Suma al array SUMA.

#open (FILE1, "<", 'coverage.txt');
open (FILE1, "<", @ARGV);#Fichero de entrada como parametro. I/O en ruta local usuario.

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  push(@suma,@fila[0]+@fila[1]);
}

close (FILE1);#Cierro la entrada del fichero fwd_rev_coverage.txt

print "\nComponiendo vector de posiciones ...\n";
#Coger las posiciones de las lineas del array SUMA que son 0 y las guarda en array POSICIONES.

my $j=0;
for (@suma) { ### |...   | [%]

  if ($suma[$j] <= $coverage_maximo)#Solo cojo aquellas lineas del fichero inferiores o iguales al coverage maximo, p.e. 0
  {
    push(@posiciones,$j+1);
  }
  $numero_bases=$j+1;
  $j++;
}                                                               

print "\nDesdoblando posiciones en iniciales y finales ...\n";
#En el array de indices POSICIONES separo los indices de inicio (a INICIO) y final (a FINAL) de cadenas-de-columna con todo "0". 

@inicio=(@posiciones[0]);#El primer valor del array INICIO tiene que ser el primero de POSICIONES y como nunca va a cumplir el if lo anado de forma manual
								
for (my $i=0;$i<=$#posiciones+1;$i++)
{	
        #Solo aquellos numeros en los que la diferencia entre ese numero y el siguiente sea mayor o igual a 2	
	if ((@posiciones[$i+1]-@posiciones[$i])>=$num_bases_entre_delecciones) 
	{	
                push(@final,@posiciones[$i]);#Guardo las coordenadas del final de la deleccion
		push(@inicio,@posiciones[$i+1]);#Guardo las coordenadas del inicio de la deleccion
	}    
}

push(@final,@posiciones[-1]);#El ultimo valor del array FINAL tiene que ser el ultimo de POSICIONES, de otra forma nunca lo obtendria
	
print "\nImprimiendo fichero de resultados ...\n";
#Imprime el principio, el final y el tamano de la deleccion en un fichero llamado output.xls (y el tamano de nuevo en el array LONGITUD).

my $tamano;#Contiene el tamano de la deleccion
open (OUT, ">", 'out_coverage.xls') or die "No puedo abrir el fichero de salida\n";

print OUT "Inicio\tFinal\tTamano\n";        

for (my $a=0;$a<=$#inicio;$a++)
{	
	$tamano=(@final[$a]-@inicio[$a])+1;#Calcula el tamaño de la deleccion
	print OUT @inicio[$a]."\t".@final[$a]."\t".(@final[$a]-@inicio[$a]+1)."\n"; 
	push (@longitud,$tamano);#meto los tamanos de las delecciones en el array
}	

print  OUT "El numero total de bases que no cambian  ", $numero_bases-sum(@longitud)," pb","\n";

close (OUT);#Cierro out_delecciones.xls

print "\nOperacion finalizada.\n\n";
