#!/usr/bin/perl
use List::Util qw (sum);

#Inputs:
#out_delecciones.xls y out_low_coverage_areas.xls
#Outputs:
#Por cada línea de ou_low_coverage_areas.xls comprobamos que no haya ningun valor de la primera columna de out_delecciones.xls en dicho rango.

my @inicios;#Array que contendra la coordenada inicial de cada deleccion.
my @inicios_excluidos;#Array que contendra los inicios excluidos del resultado.


#Creamos un vector con las posiciones inicio de out_delecciones.xls
open (FILE2, "<", 'out_delecciones.xls');

print "\nGenerando vector ...\n\n";

for (my $i=1; my $lines=<FILE2>;$i++)
{
  @fila = split (/\s+/,$lines);
    
  if ( (@fila[0] ne "Inicio") && (@fila[0] ne "El") )
  {
    push(@inicios,@fila[0]."\t".@fila[1]."\t".@fila[2]); 
  }
}

close (FILE2);

#Recorremos out_low_coverage_areas.xls y comparamos con el vector.
open (FILE1, "<", 'out_low_coverage_areas.xls');
open (OUT, ">", 'out_filtered_low_coverage_areas.xls') or die "No puedo abrir out_filtered_low_coverage_areas.xls.\n";

print "Recorriendo fichero para obtener excluidos ...\n\n";

print OUT "Inicio\tFinal\tTamano\n";

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);

  $token = 0; 
      
  for (my $j=0;$j<=$#inicios;$j++)
  {
    if ((@inicios[$j] >= @fila[0]) && (@inicios[$j] <= @fila[1]))
    {
      $token = 1;
    }
  }
  
  if (($token == 0) && ($fila[0] != "Inicio")  && ($fila[0] != "El"))
  {
    print OUT $fila[0]."\t".$fila[1]."\t".$fila[2]."\n";
  }
}

print "... Fin del Proceso\n\n";

close (FILE1);
close (OUT);#Cierro out_filtered_low_coverage_areas.xls

