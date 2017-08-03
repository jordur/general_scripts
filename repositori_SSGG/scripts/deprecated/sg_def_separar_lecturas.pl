#!/usr/bin/perl
use List::Util qw (sum);

#Sobre un fichero F3_R3.mates como este ... /data/results/Solid0065/Solid0065_20080722_1_MPControl_BACs_Wt_P2/Wt_MP_2to3Kb_005
#/results/secondary.20080811111749563/s_matching/Solid0065_20080722_1_MPControl_BACs_Wt_P2_Wt_MP_2to3Kb_005_F3.csfasta.ma.25.3
#... introducido como parametro por comando, tras desestimar las 17 lineas de cabecera, queremos generar varios ficheros:
#El primero con las lineas que no tienen coma (las que no tiene lecturas) y su sgte linea. El segundo con las que tienen 
#coma y acaban en .0, el tercero con .1, el cuarto con .2, el quinto con .3, el sexto con las que no encajen en los anteriores,
#el septimo con las que acaban en .0 y .1, el octavo con .0 .1 y .2 y el noveno con .0 .1 .2 y .3

open (FILE1, "<", @ARGV[0]);#Fichero de entrada como parametro. I/O en ruta local usuario.
open (OUTRESTO, ">", 'out_separar_lecturas_resto_lecturas.xls') or die "No puedo abrir out_separar_lecturas_resto_lecturas.xls.\n";
open (OUTSIN, ">", 'out_separar_lecturas_sin_lecturas.xls') or die "No puedo abrir out_separar_lecturas_sin_lecturas.xls.\n";
open (OUT0, ">", 'out_separar_lecturas_0_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_diferencias.xls.\n";
open (OUT1, ">", 'out_separar_lecturas_1_diferencias.xls') or die "No puedo abrir out_separar_lecturas_1_diferencias.xls.\n";
open (OUT2, ">", 'out_separar_lecturas_2_diferencias.xls') or die "No puedo abrir out_separar_lecturas_2_diferencias.xls.\n";
open (OUT3, ">", 'out_separar_lecturas_3_diferencias.xls') or die "No puedo abrir out_separar_lecturas_3_diferencias.xls.\n";
open (OUT4, ">", 'out_separar_lecturas_4_diferencias.xls') or die "No puedo abrir out_separar_lecturas_4_diferencias.xls.\n";
open (OUT5, ">", 'out_separar_lecturas_5_diferencias.xls') or die "No puedo abrir out_separar_lecturas_5_diferencias.xls.\n";
open (OUT0_1, ">", 'out_separar_lecturas_0_1_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_1_diferencias.xls.\n";
open (OUT0_1_2, ">", 'out_separar_lecturas_0_1_2_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_1_2_diferencias.xls.\n";
open (OUT0_1_2_3, ">", 'out_separar_lecturas_0_1_2_3_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_1_2_3_diferencias.xls.\n";
open (OUT0_1_2_3_4, ">", 'out_separar_lecturas_0_1_2_3_4_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_1_2_3_4_diferencias.xls.\n";
open (OUT0_1_2_3_4_5, ">", 'out_separar_lecturas_0_1_2_3_4_5_diferencias.xls') or die "No puedo abrir out_separar_lecturas_0_1_2_3_4_5_diferencias.xls.\n";
print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[0]."> ...\n";
my $destino = "";

for (my $i=1; my $lines=<FILE1>;$i++)
{	
  #if ($i > )
  {#Si linea no es cabecera
    #print "\n\n\nLa linea $i es:".$lines."\n";
    $sub = substr($lines,0,1);     
    if ($sub eq ">")
    {#Si empieza por >
      #print "SUB >\n";      
      $pos = (index($lines,","));
      if ($pos >= 0) 
      {#Si hay ,
        #print "Coma en pos-".index($lines,",")."-\n";        
        if (index($lines,",",$pos+1) >= 0)
        {#Si hay una segunda ,
          print OUTRESTO $lines;
          $destino = "RESTO";        
        }else
        {
          $punto = index($lines,".",$pos+1);
          if ($punto >= 0)
          {#Aunque siempre ha de haber un punto
            $sufijo = substr($lines,$punto+1,1);
            if ($sufijo eq "0")          
            {
              print  OUT0 $lines;
              print  OUT0_1 $lines;
              print  OUT0_1_2 $lines;
              print  OUT0_1_2_3 $lines;
	      print  OUT0_1_2_3_4 $lines;
	      print  OUT0_1_2_3_4_5 $lines;
              $destino = "0";
            }
            if ($sufijo eq "1")
            {
              print  OUT1 $lines;
              print  OUT0_1 $lines;
              print  OUT0_1_2 $lines;
              print  OUT0_1_2_3 $lines;
	      print  OUT0_1_2_3_4 $lines;
	      print  OUT0_1_2_3_4_5 $lines;
                                                                                                  
              $destino = "1";
            }
            if ($sufijo eq "2")
            {
              print  OUT2 $lines;
              print  OUT0_1_2 $lines;
              print  OUT0_1_2_3 $lines;
	      print  OUT0_1_2_3_4 $lines;
	      print  OUT0_1_2_3_4_5 $lines;
              $destino = "2";               
            }
            if ($sufijo eq "3")
            {
              print  OUT3 $lines;
              print  OUT0_1_2_3 $lines;
	      print  OUT0_1_2_3_4 $lines;
	      print  OUT0_1_2_3_4_5 $lines;
              $destino = "3";
            }
	    if ($sufijo eq "4")
            {
              print  OUT4 $lines;
	      print  OUT0_1_2_3_4 $lines;
              print  OUT0_1_2_3_4_5 $lines;
              $destino = "4";
            }
            if ($sufijo eq "5")
            {
              print  OUT5 $lines;
              print  OUT0_1_2_3_4_5 $lines;
              $destino = "5";
            }
                                                              
          }else
          {#Imposible, pero bueno ...
            $destino = "";          
          }        
        }        
      }else
      {
        print OUTSIN $lines;
        #print  $lines;
        $destino = "SIN";      
      }        
    }else #Segun lo que tengamos en la vble. destino lo metemos en un fichero, en otro o en ninguno.
    {
      #print "SUB T\n";      
      if ($destino ne "")
      {
        if ($destino eq "RESTO")
        {
          print OUTRESTO $lines;
          #print "RESTO-".$lines;
        }                                            
        if ($destino eq "SIN")
        {
          print OUTSIN $lines;
          #print "SIN-".$lines;           
        }
        if ($destino eq "0")
        {
          print OUT0 $lines;
          print OUT0_1 $lines;
          print OUT0_1_2 $lines;
          print OUT0_1_2_3 $lines;
          print OUT0_1_2_3_4 $lines;
	  print OUT0_1_2_3_4_5 $lines;
          #print "0-".$lines;
        }
        if ($destino eq "1")
        {
          print OUT1 $lines;
          print OUT0_1 $lines;
          print OUT0_1_2 $lines;
          print OUT0_1_2_3 $lines;
	  print OUT0_1_2_3_4 $lines;
	  print OUT0_1_2_3_4_5 $lines;
          #print "1-".$lines;
        }
        if ($destino eq "2")
        {
          print OUT2 $lines;
          print OUT0_1_2 $lines;
          print OUT0_1_2_3 $lines;
	  print OUT0_1_2_3_4 $lines;
	  print OUT0_1_2_3_4_5 $lines;
          #print "2-".$lines;
        }
        if ($destino eq "3")
        {
          print OUT3 $lines;
          print OUT0_1_2_3 $lines;
	  print OUT0_1_2_3_4 $lines;
	  print OUT0_1_2_3_4_5 $lines;
          #print "3-".$lines;
        }
	 if ($destino eq "4")
        {
          print OUT4 $lines;
          print OUT0_1_2_3_4 $lines;
	  print OUT0_1_2_3_4_5 $lines;
          #print "4-".$lines;
        }
	 if ($destino eq "5")
        {
          print OUT5 $lines;
          print OUT0_1_2_3_4_5 $lines;
          #print "5-".$lines;
        }

        $destino = "";        
      }#fin si hay destino
    }#fin si no empieza por >
  }#fin si linea no es cabecera
}#fin por cada linea

close (FILE1);
close (OUTRESTO);
close (OUTSIN);
close (OUT0);
close (OUT1);
close (OUT2);
close (OUT3);
close (OUT4);
close (OUT5);
close (OUT0_1);
close (OUT0_1_2);
close (OUT0_1_2_3);
close (OUT0_1_2_3_4);
close (OUT0_1_2_3_4_5);

print "\n... Fin del proceso.\n\n";
