#!/usr/bin/perl
use List::Util qw (sum);

#Inputs:
#coverage.txt (1 columna) ó fwd_rev_coverage.txt (2 columnas)
#Outputs:
#Para obtener la media del numero de colisiones hemos de calcular la suma de los valores de la columna resultante y su numero de valores.

open (FILE1, "<", @ARGV);

my $valor;#Cada fila.
$valor=0;
my $numero_valores;#Suma el numero de valores.
$numero_valores=0;
my $numero_ceros;#Sumar el numero de 0.
$numero_ceros=0;
my $sumatorio;#Suma valores.
$sumatorio=0;
my $media;#Valor medio.
$media=0;

for (my $i=1; my $lines=<FILE1>;$i++)
{
   @fila = split (/\s+/,$lines);
   $valor = $fila[0]+$fila[1];        

  if ($valor > 0)#Solo cojo aquellas lineas del fichero superiores a "0"
  {
    #print $lines,"\n";#Saco los valores a pantalla.
    $sumatorio = $sumatorio + $valor;    
  }else
  {
    $numero_ceros++;
  }
                                                                                                
  $numero_valores = $i - $numero_ceros;
}

print "Numero de valores: ".$numero_valores,"\n";#Saco el sumatorio de los valores a pantalla.
print "Sumatorio: ".$sumatorio,"\n";#Saco el numero de lineas a pantalla.
$media = $sumatorio / $numero_valores;
print "Media Aritmetica: ".$media,"\n";#Saco la media a pantalla.


close (FILE1);#Cierro la entrada del fichero coverage.txt


