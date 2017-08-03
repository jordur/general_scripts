#!/usr/bin/perl
use List::Util qw (sum);

#Inputs:
#Solid0065_20080722_1_MPControl_BACs_Wt_P2_Wt_MP_2to3Kb_005_F3_QV.qual, con los valores de calidad de las lecturas
#out_separar_lecturas_sin_lecturas.xls, o cualquiera de los ficheros out_separar_lecturas_*.xls fruto de cal_lecturas_unicas.pl
#Outputs:
#out_valores_calidad_filtrados.xls, con las parejas nombre-lectura de *.qual que estén en out_separar_lecturas_*.xls


open (FILE1, "<", @ARGV[0]);# *.qual
open (FILE2, "<", @ARGV[1]);# out_separar_lecturas_*.xls

print "\nCargando el fichero  <".$ARGV[0]."> en memoria ...\n";

for (my $i=1; my $lines=<FILE1>;$i++)
{	
  #@fila = split (/\s+/,$lines);
  push (@fichero1,$lines);
}

print "\nCargando el fichero  <".$ARGV[1]."> en memoria ...\n";

for (my $i=1; my $lines=<FILE2>;$i++)
{
  #@fila = split (/\s+/,$lines);
  #push (@fichero2,$fila[0]);
  push (@fichero2,$lines);
}

close (FILE1);
close (FILE2);

open (OUT, ">", "out_valores_calidad_filtrados.xls") or die "No puedo abrir el fichero de salida\n";

print "\nCotejando ambos ficheros ...\n";

for (my $i=0; $i<=$#fichero1;$i++)
{
  if ($token == 1)
  {
    $token = 2;
  }
  else
  {
    $token = 0;
  }
    
  if ($token == 0)
  {
    for (my $j=0;$j<=$#fichero2;$j++)
    {
      if (@fichero2[$j] eq @fichero1[$i]) 
      {
        $token = 1;
      }
    }
  }
                            
  if ($token != 0)
  {
    print OUT @fichero1[$i];
  }
}

print "\nOperacion finalizada.\n\n";
close (OUT);
