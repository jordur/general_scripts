#!/usr/bin/perl
use List::Util qw (sum);
#use Smart::Comments '###';

#Objetivo: #Determinar cuales de las lineas de out_delecciones.xls (con las zonas de la secuencia de referencia donde no mapean lecturas) son delecciones reales y cuales no. 
#Aquellas cuyos rangos se inscriban en al menos uno de los rangos de out_delecciones_f_200.xls (que son las zonas de la secuencia de referencia donde se de una alta concentracion 
#de emparejamientos con tamano del inserto mayor a 3600) seran las reales (en lugar de duplicaciones) y las pasamos a un fichero. 
#Inputs: #Un fichero tipo F3_R3.mates.non-redundant.ref1.CATEGORIA donde categoria es AAC, AAA, ... 
#Un fichero tipo out_coverage.xls (resultado de cal_coverage). 
#Diferencia Maxima por teclado, p.e. 200. 
#Posiciones Minimas por grupo, p.defecto 2. 
#Outputs: #Queremos quedarnos primero con solo las columnas F3 position y R3 position, y quitaremos el signo menos de los valores que lo tengan. Metemos estos valores en un vector 
#DESORDEN (primero el menor) junto con el tamano que resulta de restar los valores positivos de F3 y R3 (diferencia). El resultado lo ordenamos ascendentemente por la segunda 
#y luego por la primera columna y lo metemos en un vector ORDEN. A partir del vector ORDEN vamos procesando valor a valor incluyendolo en un grupo1, si la diferencia de un valor 
#con el siguiente es mayor de, p.e. 200 (esta diferencia maxima se la introducimos tb por teclado), seguimos introduciendo en grupo2, y asi sucesivamente. Teniendo en cuenta que 
#solo se contemplan los grupos formados por dos o mas posiciones (y atendiendo al valor introducido por teclado). Dentro de cada grupo establecemos nuevos grupos siguiendo el 
#mismo criterio de diferencia de valores esta vez para la segunda columna. 
#La salida ha de ser un fichero que lleva por nombre out_delecciones_f_200.xls que contenga para cada grupo; su numero de subgrupo, numero de posiciones subgrupo, 
#posiciones inicio-final del grupo contenedor, posiciones inicio-final del subgrupo, el valor medio de las posiciones de la columna primera, el valor medio de las posiciones 
#de la columna segunda y el valor medio de las diferencias (tercera columna). 
#Otro fichero out_delecciones_g_200.xls con informacion relativa a las lineas de los grupos. 
#Un tercer fichero out_delecciones_s_200.xls con informacion relativa a las lineas de los subgrupos. 
#Un cuarto fichero out_delecciones_r_200.xls con los rangos de los coverages que son delecciones reales, aunque duplicados. 
#Un quinto fichero out_delecciones_d_200.xls con los rangos de las delecciones reales tras filtrar por tamano. El resultado principal.

open (FILE1, "<", @ARGV[0]);#Fichero de entrada como parametro. I/O en ruta local usuario.
open (FILE2, "<", @ARGV[1]);#Fichero de entrada como parametro. I/O en ruta local usuario.

print "\nDiferencia Maxima: ";
$max_diferencia = <STDIN>;
chop ($max_diferencia);
print "Posiciones Minimas [2 por defecto]: ";
$min_posiciones = <STDIN>;
chop ($min_posiciones);
if ($min_posiciones <= 2)
{
  $min_posiciones = 2;
}

print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[0]."> ...\n";

my $j = 0;
for (my $i=1; my $lines=<FILE1>;$i++)
{	
  @fila = split (/\s+/,$lines);

  #Por si quiero filtrar, p.e. cojo aquellas lineas del fichero que contienen la categoria correcta
  if (1)
  {
    $valor1 = $fila[8];
    $valor2 = $fila[9];

    $posicion = index($valor1,"-");
    if ($posicion != -1)
    {
      $valor1 = substr($valor1,$posicion+1,length($valor1));  
    }

    $posicion = index($valor2,"-");
    if ($posicion != -1)
    {
      $valor2 = substr($valor2,$posicion+1,length($valor2));
    }
    
    if ($valor1 >= $valor2)
    { 
      $val [$j] = $valor2;
      $val_2 [$j] = $valor1;
      $dif [$j] = $valor1 - $valor2;
    }else
    {
      $val [$j] = $valor1;
      $val_2 [$j] = $valor2;
      $dif [$j] = $valor2 - $valor1;
    }
    
    $dupla = $val_2[$j]."|".$val[$j]."|".$dif[$j]; 
    push (@desorden,$dupla);

    $j++;
  }
}

print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[1]."> ...\n";

for (my $i=1; my $lines=<FILE2>;$i++)
{
  @fila = split (/\s+/,$lines);
  if (($fila[0] ne "Inicio") and ($fila[0] ne "El"))#Solo cojo aquellas lineas del fichero que no son cabecera o pie.
  {
    push(@out_coverage,$lines);
  }
}

#IMPRESION OUT_COVERAGE
#foreach (@out_coverage)
#{
#  print $_;
#}
    
close (FILE1);
close (FILE2);

open (G_OUT, ">", "out_delecciones_g_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";
open (S_OUT, ">", "out_delecciones_s_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";
open (F_OUT, ">", "out_delecciones_f_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";
open (R_OUT, ">", "out_delecciones_r_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";
open (D_OUT, ">", "out_delecciones_d_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";

print "\nOrdenando la informacion seleccionada ...\n";

@orden = sort { $a <=> $b } @desorden;

#IMPRESION ORDEN CRITERIO 2
foreach (@orden)
{
  #print $_."\n";
  ($col_2,$col_1,$dif_3)=split(/\|/,$_);
  $dupla = $col_1."|".$col_2."|".$dif_3;      
  push (@desorden_2,$dupla);
}

@orden = sort { $a <=> $b } @desorden_2;

#IMPRESION ORDEN CRITERIO 1
#open (OUT2, ">", "out_pantalla_".$max_diferencia.".xls") or die "No puedo abrir el fichero de salida\n";
#print OUT2 "Col_1\tCol_2\tDiferencia\n";
#print "Col_1\tCol_2\tDiferencia\n";
#foreach (@orden)
#{
#  print OUT2 $_."\n";
#}
#close (OUT2);

#MATRIZ BIDIMENSIONAL EXTINTA
#foreach my $k (0 .. $j-1)
#{
#  foreach my $l (0 .. 1)
#  {
#    print " Madre k: $k l: $l $madre[$k][$l]";
#  }
#  print "\n";
#}

print "\nNotacion:\n";
print "g - grupos\n";
print "s - subgrupos\n";
print "f - final\n";
print "r - out_coverage reales\n";
print "d - delecciones reales\n";

print "\nExportando los resultados a los ficheros <out_delecciones_[g-s-f-r-d]_".$max_diferencia.".xls> ...\n";

print G_OUT "Tipo_Grupo\tNum_Grupo\tNum_Posiciones\tRango_1\n";
print S_OUT "Tipo_Grupo\tNum_Grupo\tNum_Posiciones\tRango_1\n";
#print F_OUT "Num_Grupo\tNum_Posiciones\tRango_1\tRango_2\tMedia_Posiciones_1\tMedia_Posiciones_2\tDistancia_Media\n";
print F_OUT "Num_Posiciones\tRango_1\tRango_2\tDistancia_Entre_Rangos\n";
print R_OUT "Inicio\tFinal\tTamano\n";
print D_OUT "Num_Posiciones\tRango_1\tRango_2\tDistancia_Entre_Rangos\n";

$num_duplas = 1;#Como maximo sera de j.
$k = 0;#Numero de vectores resultantes.
$l = 0;#Numero de valores del vector en curso.
$posicion = 0;
$posicion_2 = 0;
$diferencia = 0;
$posicion_previa = 0;
$num_grupo = 0;

foreach $dupla(@orden) {  ### |...    | [%]

  $l++;#Numero de valores del vector en curso.
  ($posicion,$posicion_2,$diferencia)=split(/\|/,$dupla);
  if ($num_duplas == 1)
  {
    $posicion_inicio = $posicion;
    $posicion_previa = $posicion;
    my @g1 = ($posicion);
    my @g2 = ($posicion_2);
    my @g3 = ($diferencia);
  }
  
  if ( ($num_duplas != 1) && (($num_duplas == $j) || ($posicion - $posicion_previa >= $max_diferencia)) )
  {
    $k++;
    
    if ($posicion - $posicion_previa >= $max_diferencia)
    {
      if ($l-1 >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
      {
        print G_OUT "GRUPO_A_PRIORI_1 ".$k."\t".($l-1)."\t".$posicion_inicio."-".$posicion_previa."\n";
      
      
# ----------- PROCESADO GRUPO BLOQUE 1
      

        #IMPRESION GRUPO
        #print "\nCCol_1\tCol_2\tDiferencia\n";
        $vaciar = @s_desorden;
        for (my $n=1; $n<=$vaciar; $n++)
        {
          $vacio = shift(@s_desorden);
        }
        foreach my $m (0 .. $#g1)
        {
          if ($l-1 >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
          {
            print G_OUT $g1[$m]."\t".$g2[$m]."\t".$g3[$m]."\n";
          }
          $s_dupla = $g2[$m]."|".$g1[$m]."|".$g3[$m];              
          push (@s_desorden,$s_dupla);              
        }      
        $vaciar = @s_orden;
        for (my $n=1; $n<=$vaciar; $n++)
        {
          $vacio = shift(@s_orden);
        } 
        @s_orden = sort { $a <=> $b } @s_desorden;
                
        #IMPRESION ORDEN CRITERIO 2
        $s_num_duplas = 1;#Como maximo sera de @s_orden.
        $s_k = 0;#Numero de vectores resultantes.
        $s_posicion = 0;
        $s_posicion_2 = 0;
        $s_diferencia = 0;
        $s_posicion_previa = 0;
        $s_j = @s_orden;
  
        #print "\nCCol_2\tCol_1\tDiferencia\n";            
        foreach $s_dupla(@s_orden)
        {#inicio for
          #print G_OUT $s_dupla."\n";
  
          ($s_posicion,$s_posicion_2,$s_diferencia)=split(/\|/,$s_dupla);
          if ($s_num_duplas == 1)
          {
            $s_posicion_inicio = $s_posicion;
            $s_posicion_previa = $s_posicion;
            my @s_g1 = ($s_posicion);
            my @s_g2 = ($s_posicion_2);
            my @s_g3 = ($s_diferencia);
          }
  
          if ( ($s_num_duplas != 1) && (($s_num_duplas == $s_j) || ($s_posicion - $s_posicion_previa >= $max_diferencia)) )
          {
            $s_k++;
  
            if ($s_posicion - $s_posicion_previa >= $max_diferencia)
            {
              $num_grupo++;
              my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
              my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
              my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
              $s_vaciar = @s_g1;
              #print S_OUT $s_vaciar."\n";
              
              @s_g2_ordenado = sort { $a <=> $b } @s_g2;
              
              if ($s_vaciar >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
              {
                print S_OUT "SUBGRUPO_1 ".$num_grupo."\t".$s_vaciar."\t".$s_posicion_inicio."-".$s_posicion_previa."\n";
                #print F_OUT $num_grupo."\t".$s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($suma_s_g2/$s_vaciar)."\t".int($suma_s_g1/$s_vaciar)."\t".int($suma_s_g3/$s_vaciar)."\n";
                print F_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
              
                my $suma_coverage_tamano = 0;
                for (my $p=1; $p<=@out_coverage; $p++)
                {
                  $out_coverage_fila = $out_coverage[$p];
                  ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                  if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion_previa) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion_previa) )
                  {
                    print F_OUT $out_coverage_fila;
                    print R_OUT $out_coverage_fila;
                    #splice(@out_coverage,$p,1); redo;
                    $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                  }
                }
                
                if ($suma_coverage_tamano > int(($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                {
                  print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
                }
              } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
                         
              foreach my $o (0..($s_vaciar-1))
              {
                print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
              }
              for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
              {
                $s_vacio = shift(@s_g1);
                $s_vacio = shift(@s_g2);
                $s_vacio = shift(@s_g3);
              }
              push (@s_g1,$s_posicion);
              push (@s_g2,$s_posicion_2);
              push (@s_g3,$s_diferencia);
  
              if ($s_num_duplas == $s_j)
              {# Como es un grupo de una unica posicion al final se descartara, pero bueno ...
                $num_grupo++;
                my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
                my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
                my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
                $s_vaciar = @s_g1;
                              
                @s_g2_ordenado = sort { $a <=> $b } @s_g2;
  
                if (0) #Solo se contemplan los grupos formados por mas de dos posiciones.
                {
                  print S_OUT "SUBGRUPO_1B ".$num_grupo."\t1\t".$s_posicion."-".$s_posicion."\n";
                  #print F_OUT $num_grupo."\t1\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[0]."\t".$s_posicion."-".$s_posicion."\t".int($suma_s_g2)."\t".int($suma_s_g1)."\t".int($suma_s_g3)."\n";            
                  print F_OUT "1\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[0]."\t".$s_posicion."-".$s_posicion."\t".int($s_posicion-$s_g2_ordenado[0])."\n";
  
                  my $suma_coverage_tamano = 0;
                  for (my $p=1; $p<=@out_coverage; $p++)
                  {
                    $out_coverage_fila = $out_coverage[$p];
                    ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                    if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion) )
                    {
                      print F_OUT $out_coverage_fila;
                      print R_OUT $out_coverage_fila;
                      #splice(@out_coverage,$p,1); redo;
                      $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                    }
                  }
                
                  if ($suma_coverage_tamano > int(($s_posicion-$s_g2_ordenado[0])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                  {
                    print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion-$s_g2_ordenado[0])."\n";
                  }
                } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
                
                foreach my $o (0..($s_vaciar-1))
                {
                  print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
                }
                for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
                {
                  $s_vacio = shift(@s_g1);
                  $s_vacio = shift(@s_g2);
                  $s_vacio = shift(@s_g3);
                }
              }
            }else
            {
              if ($s_num_duplas == $s_j)
              {
                $num_grupo++;
  
                push (@s_g1,$s_posicion);
                push (@s_g2,$s_posicion_2);
                push (@s_g3,$s_diferencia);
                                
                my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
                my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
                my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
                $s_vaciar = @s_g1;
                #print S_OUT $s_vaciar."\n";
  
                @s_g2_ordenado = sort { $a <=> $b } @s_g2;
  
                if ($s_vaciar >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
                {                                          
                  print S_OUT "SUBGRUPO_2 ".$num_grupo."\t".$s_vaciar."\t".$s_posicion_inicio."-".$s_posicion."\n";
                  #print F_OUT $num_grupo."\t".$s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion."\t".int($suma_s_g2/$s_vaciar)."\t".int($suma_s_g1/$s_vaciar)."\t".int($suma_s_g3/$s_vaciar)."\n";
                  print F_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
  
                  my $suma_coverage_tamano = 0;
                  for (my $p=1; $p<=@out_coverage; $p++)
                  {
                    $out_coverage_fila = $out_coverage[$p];
                    ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                    if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion) )
                    {
                      print F_OUT $out_coverage_fila;
                      print R_OUT $out_coverage_fila;
                      #splice(@out_coverage,$p,1); redo;
                      $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                    }
                  }
                
                  if ($suma_coverage_tamano > int(($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                  {
                    print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
                  }
                } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
                
                foreach my $o (0..($s_vaciar-1))
                {
                  print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
                }
                for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
                {
                  $s_vacio = shift(@s_g1);
                  $s_vacio = shift(@s_g2);
                  $s_vacio = shift(@s_g3);
                }
              }
            }
  
            $s_posicion_inicio = $s_posicion;
          }else
          {
            push (@s_g1,$s_posicion);
            push (@s_g2,$s_posicion_2);
            push (@s_g3,$s_diferencia);
          }
  
          $s_num_duplas++;
          $s_posicion_previa = $s_posicion;
                                
        }#fin for
                      
                            
# ----------- FIN PROCESADO GRUPO BLOQUE 1
      
      }#fin Solo se contemplan los grupos formados por mas de dos posiciones.
      
      $vaciar = @g1;
      for (my $n=1; $n<=$vaciar; $n++)
      {
        $vacio = shift(@g1);
        $vacio = shift(@g2);
        $vacio = shift(@g3);
      }
      my @g1 = ($posicion);
      my @g2 = ($posicion_2);
      my @g3 = ($diferencia);
      
      if ($num_duplas == $j)
      {# Como es un grupo de una unica posicion al final se descartara, pero bueno ...        
        if (0) #Solo se contemplan los grupos formados por mas de dos posiciones.
        {        
          print G_OUT "GRUPO_A_PRIORI_1B ".($k+1)."\t1\t".$posicion."-".$posicion."\n"; 

          # ----------- PROCESADO GRUPO BLOQUE 1B

          #IMPRESION GRUPO
          #print "\nCol_1\tCol_2\tDiferencia\n";
          foreach my $m (0 .. $#g1)
          {
            print G_OUT $g1[$m]."\t".$g2[$m]."\t".$g3[$m]."\n";
          }

          # ----------- FIN PROCESADO GRUPO BLOQUE 1B

        } #fin Solo se contemplan los grupos formados por mas de dos posiciones.        
      }
      
    }else
    {
      if ($num_duplas == $j)
      {
        if ($l >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
        {
          print G_OUT "GRUPO_A_PRIORI_2 ".$k."\t".$l."\t".$posicion_inicio."-".$posicion."\n";


# ----------- PROCESADO GRUPO BLOQUE 2 ---------- ATENCION COPY-PASTE 


          #IMPRESION GRUPO
          #print "\nCCol_1\tCol_2\tDiferencia\n";
          $vaciar = @s_desorden;
          for (my $n=1; $n<=$vaciar; $n++)
          {
            $vacio = shift(@s_desorden);
          }
          foreach my $m (0 .. $#g1)
          {
            print G_OUT $g1[$m]."\t".$g2[$m]."\t".$g3[$m]."\n";
            $s_dupla = $g2[$m]."|".$g1[$m]."|".$g3[$m];              
            push (@s_desorden,$s_dupla);              
          }  
          print G_OUT $posicion."\t".$posicion_2."\t".$diferencia."\n";
          $vaciar = @s_orden;
          for (my $n=1; $n<=$vaciar; $n++)
          {
            $vacio = shift(@s_orden);
          } 
          @s_orden = sort { $a <=> $b } @s_desorden;
                  
          #IMPRESION ORDEN CRITERIO 2
          $s_num_duplas = 1;#Como maximo sera de @s_orden.
          $s_k = 0;#Numero de vectores resultantes.
          $s_posicion = 0;
          $s_posicion_2 = 0;
          $s_diferencia = 0;
          $s_posicion_previa = 0;
          $s_j = @s_orden;
    
          #print "\nCCol_2\tCol_1\tDiferencia\n";            
          foreach $s_dupla(@s_orden)
          {#inicio for
            #print G_OUT $s_dupla."\n";
    
            ($s_posicion,$s_posicion_2,$s_diferencia)=split(/\|/,$s_dupla);
            if ($s_num_duplas == 1)
            {
              $s_posicion_inicio = $s_posicion;
              $s_posicion_previa = $s_posicion;
              my @s_g1 = ($s_posicion);
              my @s_g2 = ($s_posicion_2);
              my @s_g3 = ($s_diferencia);
            }
    
            if ( ($s_num_duplas != 1) && (($s_num_duplas == $s_j) || ($s_posicion - $s_posicion_previa >= $max_diferencia)) )
            {
              $s_k++;
    
              if ($s_posicion - $s_posicion_previa >= $max_diferencia)
              {
                $num_grupo++;
                my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
                my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
                my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
                $s_vaciar = @s_g1;
                #print S_OUT $s_vaciar."\n";
                
                @s_g2_ordenado = sort { $a <=> $b } @s_g2;
  
                if ($s_vaciar >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
                {
                  print S_OUT "SUBGRUPO_1 ".$num_grupo."\t".$s_vaciar."\t".$s_posicion_inicio."-".$s_posicion_previa."\n";
                  #print F_OUT $num_grupo."\t".$s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($suma_s_g2/$s_vaciar)."\t".int($suma_s_g1/$s_vaciar)."\t".int($suma_s_g3/$s_vaciar)."\n";
                  print F_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
  
                  my $suma_coverage_tamano = 0;
                  for (my $p=1; $p<=@out_coverage; $p++)
                  {
                    $out_coverage_fila = $out_coverage[$p];
                    ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                    if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion_previa) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion_previa) )
                    {
                      print F_OUT $out_coverage_fila;
                      print R_OUT $out_coverage_fila;                                        
                      #splice(@out_coverage,$p,1); redo;
                      $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                    }
                  }
                
                  if ($suma_coverage_tamano > int(($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                  {
                    print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
                  }
                } #fin Solo se contemplan los grupos formados por mas de dos posiciones.        
        
                foreach my $o (0..($s_vaciar-1))
                {
                  print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
                }
                for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
                {
                  $s_vacio = shift(@s_g1);
                  $s_vacio = shift(@s_g2);
                  $s_vacio = shift(@s_g3);
                }
                push (@s_g1,$s_posicion);
                push (@s_g2,$s_posicion_2);
                push (@s_g3,$s_diferencia);
    
                if ($s_num_duplas == $s_j)
                {# Como es un grupo de una unica posicion al final se descartara, pero bueno ...
                  $num_grupo++;
                  my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
                  my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
                  my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
                  $s_vaciar = @s_g1;
    
                  @s_g2_ordenado = sort { $a <=> $b } @s_g2;
  
                  if (0) #Solo se contemplan los grupos formados por mas de dos posiciones.
                  {
                    print S_OUT "SUBGRUPO_1B ".$num_grupo."\t1\t".$s_posicion."-".$s_posicion."\n";
                    #print F_OUT $num_grupo."\t1\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[0]."\t".$s_posicion."-".$s_posicion."\t".int($suma_s_g2)."\t".int($suma_s_g1)."\t".int($suma_s_g3)."\n";
                    print F_OUT "1\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[0]."\t".$s_posicion."-".$s_posicion."\t".int($s_posicion-$s_g2_ordenado[0])."\n";            
  
                    my $suma_coverage_tamano = 0;
                    for (my $p=1; $p<=@out_coverage; $p++)
                    {
                      $out_coverage_fila = $out_coverage[$p];
                      ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                      if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion) )
                      {
                        print F_OUT $out_coverage_fila;
                        print R_OUT $out_coverage_fila;                                          
                        #splice(@out_coverage,$p,1); redo;
                        $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                      }
                    }
                
                    if ($suma_coverage_tamano > int(($s_posicion-$s_g2_ordenado[0])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                    {
                      print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion-$s_g2_ordenado[0])."\n";
                    }
                  } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
    
                  foreach my $o (0..($s_vaciar-1))
                  {
                    print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
                  }
                  for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
                  {
                    $s_vacio = shift(@s_g1);
                    $s_vacio = shift(@s_g2);
                    $s_vacio = shift(@s_g3);
                  }
                }
              }else
              {
                if ($s_num_duplas == $s_j)
                {
                  $num_grupo++;
    
                  push (@s_g1,$s_posicion);
                  push (@s_g2,$s_posicion_2);
                  push (@s_g3,$s_diferencia);
                                  
                  my $suma_s_g1; map {$suma_s_g1 += $_} @s_g1;
                  my $suma_s_g2; map {$suma_s_g2 += $_} @s_g2;
                  my $suma_s_g3; map {$suma_s_g3 += $_} @s_g3;
  
                  $s_vaciar = @s_g1;
                  #print S_OUT $s_vaciar."\n";
    
                  @s_g2_ordenado = sort { $a <=> $b } @s_g2;
  
                  if ($s_vaciar >= $min_posiciones) #Solo se contemplan los grupos formados por mas de dos posiciones.
                  {
                    print S_OUT "SUBGRUPO_2 ".$num_grupo."\t".$s_vaciar."\t".$s_posicion_inicio."-".$s_posicion."\n";
                    #print F_OUT $num_grupo."\t".$s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion."\t".int($suma_s_g2/$s_vaciar)."\t".int($suma_s_g1/$s_vaciar)."\t".int($suma_s_g3/$s_vaciar)."\n";
                    print F_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
  
                    my $suma_coverage_tamano = 0;
                    for (my $p=1; $p<=@out_coverage; $p++)
                    {
                      $out_coverage_fila = $out_coverage[$p];
                      ($out_coverage_inicio,$out_coverage_final,$out_coverage_tamano)=split(/\s+/,$out_coverage_fila);
                      if ( ($out_coverage_inicio >= $s_g2_ordenado[0]) and ($out_coverage_inicio <= $s_posicion) and ($out_coverage_final >= $s_g2_ordenado[0]) and ($out_coverage_final <= $s_posicion) )
                      {
                        print F_OUT $out_coverage_fila;
                        print R_OUT $out_coverage_fila;
                        #splice(@out_coverage,$p,1); redo;
                        $suma_coverage_tamano = $suma_coverage_tamano + $out_coverage_tamano;
                      }
                    }
                
                    if ($suma_coverage_tamano > int(($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])/3)) #Solo se contemplan los grupos donde la suma de los tamanos sea mayor del tercio de la diferencia entre rangos.
                    {
                      print D_OUT $s_vaciar."\t".$s_g2_ordenado[0]."-".$s_g2_ordenado[@s_g2_ordenado-1]."\t".$s_posicion_inicio."-".$s_posicion_previa."\t".int($s_posicion_inicio-$s_g2_ordenado[@s_g2_ordenado-1])."\n";
                    }
                  } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
                            
                  foreach my $o (0..($s_vaciar-1))
                  {
                    print S_OUT @s_g2[$o]."\t".@s_g1[$o]."\t".@s_g3[$o]."\n";
                  }
                  for (my $s_n=1; $s_n<=$s_vaciar; $s_n++)
                  {
                    $s_vacio = shift(@s_g1);
                    $s_vacio = shift(@s_g2);
                    $s_vacio = shift(@s_g3);
                  }
                }
              }
    
              $s_posicion_inicio = $s_posicion;
            }else
            {
              push (@s_g1,$s_posicion);
              push (@s_g2,$s_posicion_2);
              push (@s_g3,$s_diferencia);
            }
    
            $s_num_duplas++;
            $s_posicion_previa = $s_posicion;
                                  
          }#fin for


# ----------- FIN PROCESADO GRUPO BLOQUE 2 ---------- ATENCION COPY-PASTE


        } #fin Solo se contemplan los grupos formados por mas de dos posiciones.
        
        $vaciar = @g1;
        for (my $n=1; $n<=$vaciar; $n++)
        {
          $vacio = shift(@g1);
          $vacio = shift(@g2);
          $vacio = shift(@g3);
        }
        my @g1 = ($posicion);
        my @g2 = ($posicion_2);
        my @g3 = ($diferencia);
      }
    }

    $l = 1;#Numero de valores del vector en curso.
    $posicion_inicio = $posicion;
    push (@g1,$posicion);
    push (@g2,$posicion_2);
    push (@g3,$diferencia);
  }else
  {
    push (@g1,$posicion);
    push (@g2,$posicion_2);
    push (@g3,$diferencia);
  }
 
  $num_duplas++;    
  $posicion_previa = $posicion;
}

print "\nOperacion finalizada.\n\n";
close (G_OUT);
close (S_OUT);
close (F_OUT);
close (R_OUT);
