#!/usr/bin/perl
use List::Util qw (sum);
use strict;
#Inputs:
#Fichero con la primera columna numerica
#Outputs:
#Resultado en pantalla del sumatorio de los valores de esa columna

open (FILE1, "<", @ARGV);

my $valor = 0;
my $sumatorio = 0;
my @fila;

for (my $i=1; my $lines=<FILE1>;$i++)
{
   @fila = split (/\s+/,$lines);
   $valor = $fila[0];        
   $sumatorio = $sumatorio + $valor;    
}

print $sumatorio,"\n";#Saco el numero de lineas a pantalla. 

close (FILE1);#Cierro la entrada del fichero
