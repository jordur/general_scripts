#!/usr/bin/perl
use List::Util qw (sum);
use strict;

#Inputs:
#Fichero con varias secuencias de referencia tipo:
#>hola
#ATAATTTTGCAGTTGACAAACCATTGAGGCTTTATAAGAGGCTCCACTAC
#ATCATTGCTTCTTGAACAAATCCCAAGGCGCATCTCATTATTCTTAGCAC
#CCAGTTTCTGATTATCCTCATCCACAACAGCATCCAAAGCCTGCTTGGCC
#>mundo
#CCAGTTTCTGATTATCCTCATCCACAACAGCATCCAAAGCCTGCTTGGCC
#ATAATTTTGCAGTTGACAAACCATTGAGGCTTTATAAGAGGCTCCACTAC
#...
#Otro fichero tipo snps.txt con los snps de la secuencia_union
#Outputs:
#out_secuencia_union.xls con la concatenacion de todas las anterioes.
#out_fasta.xls con nombre secuencia de referencia, posici inicio, final y tamano.
#ID_ref	Start     End      Size
#hola	1       50	50
#mundo	51	120     70
#Solo si existe el fichero snps.txt de entrada, out_snps_union.xls con el formato:
#ID_ref    # cov    ref     consen     score    confi     Single   Paired   score    conf      Single   Paired   coord_in_ref      corrd_in_seq
#hola	135       G	T       0.0000  0.0000  0/0	0/0     1.0000  0.9444  57/31    78/45	2       2                        
#hola   219       G     A       0.0046  0.9425  1/1     0/0     0.9774  0.9333  89/48   125/60  45      45
#mundo  198       T     C	0.0000  0.0000  0/0     0/0	0.9951  0.9479  86/44   111/57  58	8


#Extraemos info fragmentos.

print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[0]."> ...\n";

open (FILE1, "<", @ARGV[0]);#Fichero de entrada como parametro. I/O en ruta local usuario.
open (OUT_S, ">", 'out_secuencia_union.xls') or die "No puedo abrir el fichero.\n";
open (OUT_F, ">", 'out_fasta.xls') or die "No puedo abrir el fichero.\n";

my @fichero1;

for (my $i=1; my $lines=<FILE1>; $i++)
{
  chomp($lines);
  push (@fichero1,$lines);
}

close (FILE1);#Cierro la entrada del fichero

print "\nProcesando las secuencias ...\n";

my @fila;
my $nom_sec;
my $secuencia_union_linea;
my @secuencia_union;
my $fasta;
my $posicion1=0;
my $posicion2=0;
my $anterior=0;

print OUT_F "ID_ref\tStart\tEnd\tSize\n";

for (my $i=0; $i<=$#fichero1; $i++)
{
  @fila = split (/\s+/,$fichero1[$i]);
  if (substr($fila[0],(0),1) eq ">")
  {      
    if ($i!=0)
    {
      $secuencia_union_linea = $nom_sec."\t".$posicion1."\t".$posicion2."\t".($posicion2-$posicion1+1);
      push (@secuencia_union,$secuencia_union_linea);
    }

    $posicion1 = ($posicion2+1);
    $nom_sec = substr($fila[0],1);
    $secuencia_union_linea = "";
    $anterior = 1;                
  }else
  {
    $fasta = $fasta.$fila[0];
    $posicion2 = $posicion2 + length($fila[0]);

    if ($i==$#fichero1)
    {
      $secuencia_union_linea = $nom_sec."\t".$posicion1."\t".$posicion2."\t".($posicion2-$posicion1+1);
      push (@secuencia_union,$secuencia_union_linea);
    }
  }
}

print OUT_S $fasta."\n";


foreach (@secuencia_union)
{
  print OUT_F $_."\n";
}

close (OUT_S);#Cierro la salida a fichero
close (OUT_F);#Cierro la salida a fichero


#Extraemos info snps.

print "\nExtrayendo la informacion seleccionada del fichero <".$ARGV[1]."> (si existe) ...\n";

open (FILE2, "<", @ARGV[1]);#Fichero de entrada como parametro. I/O en ruta local usuario.
open (OUT_U, ">", 'out_snps_union.xls') or die "No puedo abrir el fichero.\n";

my $linea_snps;
my @snps;
my $valor1;
my $valor2;
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
my $c10;
my $c11;
my $c12;
my $valor1;
my $posicion1;
my $posicion2;

for (my $i=1; my $lines=<FILE2>;$i++)
{
  @fila = split (/\s+/,$lines);
  ($c0,$c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9,$c10,$c11)=split(/\s+/,$lines);

  if (($c0 ne "#") && ($c0 ne ""))
  {
    foreach (@secuencia_union)
    {
      my $e0;
      my $e1;
      my $e2;
      my $e3;
      ($e0,$e1,$e2,$e3)=split(/\t/,$_);
    
      if (($e1<=$c11) and ($e2>=$c11))
      {
        $valor1 = $e0;
        $valor2 = $c11-$e1+1;
        #print $e0." ".$e1." ".$e2." ".$e3." c11-".$c11."\n";
      }
    }
  
    $linea_snps = $valor1."\t".$c0."\t".$c1."\t".$c2."\t".$c3."\t".$c4."\t".$c5."\t".$c6."\t".$c7."\t".$c8."\t".$c9."\t".$c10."\t".$c11."\t".$valor2;
    push (@snps,$linea_snps);
  }
}   

print OUT_U "#ID_ref\tcov\tref\tconsen\tscore\tconfi\tSingle\tPaired\tscore\tconf\tSingle\tPaired\tcoord_in_ref\tcoord_in_seq\n";

foreach (@snps)
{
  print OUT_U $_."\n";
}

print "\nOperacion finalizada.\n\n";

close (OUT_U);#Cierro la salida a fichero
close (FILE2);#Cierro la entrada del fichero
