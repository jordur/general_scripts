#!/usr/bin/perl
use List::Util qw (sum);
use strict;

#Inputs:
#counts.txt
#Outputs:
#out_wt_coverage.xls

my @fila;
my $linea2;
my $linea3;
my $linea4;
my @desorden;
my @orden;
my @orden2;
my $divisor;

print "\nDivisor [1 por defecto]: ";
$divisor = <STDIN>;
chop ($divisor);

if ($divisor <= 1)
{
  $divisor = 1;  
}

#Metemos cada linea en el array, menos las que tengan 0 en col 7.

print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[0]."> ...\n";

#open (FILE1, "<", 'counts.txt');
open (FILE1, "<", @ARGV);#Fichero de entrada como parametro. I/O en ruta local usuario.

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  if ($fila[7]!=0)
  {
    $linea2 = $fila[8]."|".$fila[0]."|".$fila[1]."|".$fila[2]."|".$fila[3]."|".$fila[4]."|".$fila[5]."|".$fila[6]."|".$fila[7]; 
    push (@desorden,$linea2);
  }
}

close (FILE1);#Cierro la entrada del fichero

#Ordenamos por ID creciente incluyendo columna ($7*35)/($4-$3)/divisor 

@orden = sort(@desorden);

my $contador = 0;
my $id_anterior;
my $id_actual;
my $suma_grupo = 0;
my $contador_grupo = 0;
my $media_grupo = 0;
my $suma_total = 0;
my $contador_total = 0;
my $media_total = 0;
my $n_elementos = @orden;

#Calculamos la media de cada grupo y la media total

foreach (@orden)
{
  $contador++;
  my $c0;
  my $c1;
  my $c2;
  my $c3;
  my $c4;
  my $c5;
  my $c6;
  my $c7;
  my $c8;
  my $valor1;
  my $posicion1;
  my $posicion2;
  ($c8,$c0,$c1,$c2,$c3,$c4,$c5,$c6,$c7)=split(/\|/,$_);
  $valor1 = sprintf("%.2f",($c7*35)/($c4-$c3)/$divisor);
  $linea3 = $c0."|".$c1."|".$c2."|".$c3."|".$c4."|".$c5."|".$c6."|".$c7."|".$c8."|".$valor1;      
  
  if ($contador == 1)
  {
    $posicion1 = index($c8,"=");
    $posicion2 = index($c8,";");
    $id_anterior = substr($c8,$posicion1+1,$posicion2-$posicion1-1);

    $suma_grupo = $valor1;
    $contador_grupo = 1;
  
    push (@orden2,$linea3);
  }else
  {#resto de lineas
  
    $posicion1 = index($c8,"=");
    $posicion2 = index($c8,";");
    $id_actual = substr($c8,$posicion1+1,$posicion2-$posicion1-1);

    if ($contador != $n_elementos)
    {

      if ($id_anterior ne $id_actual)
      {#cambia id
        $media_grupo = sprintf("%.2f",$suma_grupo/$contador_grupo);
        push (@orden2,"Media_Grupo:|".$media_grupo);

        $id_anterior = $id_actual;
        $suma_grupo = $valor1;
        $contador_grupo = 1;    
        $suma_total = $suma_total + $media_grupo;
        $contador_total++;
      }else
      {#no cambia id
        $suma_grupo = $suma_grupo + $valor1;
        $contador_grupo++;  
      }#fin no cambia id
  
      push (@orden2,$linea3);
  
    }else
    {

      if ($id_anterior ne $id_actual)
      {#cambia id
        $media_grupo = sprintf("%.2f",$suma_grupo/$contador_grupo);
        push (@orden2,"Media_Grupo:|".$media_grupo);
        push (@orden2,$linea3);
        push (@orden2,"Media_Grupo:|".sprintf("%.2f",$valor1));

        $suma_total = $suma_total + $media_grupo;
        $contador_total++;
      }else
      {#no cambia id
        push (@orden2,$linea3);
        
        $suma_grupo = $suma_grupo + $valor1;
        $contador_grupo++;

        $media_grupo = sprintf("%.2f",$suma_grupo/$contador_grupo);
        push (@orden2,"Media_Grupo:|".$media_grupo);
                
      }#fin no cambia id
    
    }
    
  }#fin resto de lineas
  
}

#print "Numero de valores: ".$contador_total,"\n";
#print "Sumatorio: ".$suma_total,"\n";
$media_total = sprintf("%.2f",$suma_total/$contador_total);
#print "Media Aritmetica: ".$media_total,"\n";

push (@orden2,"Media_Total:|".$media_total);

#Imprime resultados y calcula la desviacion tipica      

open (OUT, ">", 'out_wt_coverage.xls') or die "No puedo abrir el fichero.\n";

my $desviacion_tipica = 0;

foreach (@orden2)
{
  my $c0;
  my $c1;
  my $c2;
  my $c3;
  my $c4;
  my $c5;
  my $c6;
  my $c7;
  my $c8;
  my $c9;
  ($c0,$c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9)=split(/\|/,$_);

  if ($c0 eq "Media_Grupo:")
  {
    $desviacion_tipica = $desviacion_tipica + ($c1-$media_total)**2;
  }
  
  $linea4 = $c0."\t".$c1."\t".$c2."\t".$c3."\t".$c4."\t".$c5."\t".$c6."\t".$c7."\t".$c8."\t".$c9;
  print OUT $linea4."\n";
}

$desviacion_tipica = sprintf("%.2f",sqrt($desviacion_tipica/$contador_total));
print OUT "Desviacion_Standard:\t".$desviacion_tipica."\n\n";#Saco la desviacion tipica.

print "\nOperacion finalizada.\n\n";

close (OUT);#Cierro la salida a fichero
