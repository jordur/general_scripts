#!/usr/bin/perl
use List::Util qw (sum);

#open (FILE1, "<", 'coverage.txt');
open (FILE1, "<", @ARGV);

#Inputs:
#coverage.txt (1 columna) ó fwd_rev_coverage.txt (2 columnas)
#Outputs:
#Ha de coger el coverage minimo y el coverage maximo para que solo coja las posiciones que tengan un coverage que este entre estos dos valores.
#Para obtener la media del numero de colisiones hemos de calcular la suma de los valores de la columna resultante y su numero de valores.
#Para obtener la desviacion tipica aplicamos la formula dt=sqrt((sumatorio(xi-media)^2)/N), donde xi es cada uno de los valores,
#y N es el numero total de valores tomados.

print "\nCoverage minimo:","\n";
$coverage_minimo=<STDIN>;
print "Media esperada:","\n";
$media_esperada=<STDIN>;
print "Coverage maximo:","\n";
$coverage_maximo=<STDIN>;

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
my $desviacion_tipica;
$desviacion_tipica=0;

print "\nCalculando media y desviacion standard ...","\n\n";

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  $valor = $fila[0]+$fila[1];        

  if (($valor > $coverage_minimo) and ($valor < $coverage_maximo))#Solo cojo aquellas lineas del fichero entre el valor minimo (=5) y el máximo (=triple de la media esperada)
  {
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

open (FILE1, "<", @ARGV);
for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  $valor = $fila[0]+$fila[1];

  if (($valor > $coverage_minimo) and ($valor < $coverage_maximo))#Solo cojo aquellas lineas del fichero entre el valor minimo (=5) y el máximo (=triple de la media esperada)
  {
    $desviacion_tipica = $desviacion_tipica + ($valor-$media_esperada)**2;
  }
}

$desviacion_tipica = sqrt($desviacion_tipica/$numero_valores);
print "Desviacion Standard: ".$desviacion_tipica,"\n\n";#Saco la desviacion tipica.
close (FILE1);#Cierro la entrada del fichero coverage.txt


