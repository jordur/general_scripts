#!/usr/bin/perl
use List::Util qw (sum);

#Inputs:
#snps.txt
#Outputs:
#Hay que sumar los numeradores entre si, y por otro lado los denominadores entre si, de las supuestas dos ultimas 
#columnas-fracciones del fichero prueba_snps.txt

my @fila;#Del prueba_snps.txt
my @primera;#Una fraccion.
my @segunda;#La otra fraccion.
my @suma1;# De los numeradores.
my @suma2;# De los denominadores.
my $posicion1;
my $posicion2;


#Metemos cada fraccion en un array.

#open (FILE1, "<", 'snps.txt');
open (FILE1, "<", @ARGV);#Fichero de entrada como parametro. I/O en ruta local usuario.

for (my $i=1; my $lines=<FILE1>;$i++)
{
  @fila = split (/\s+/,$lines);
  push(@primera,@fila[9]);
  push(@segunda,@fila[10]);  
}

close (FILE1);#Cierro la entrada del fichero


#Coger los numeradores y sumarlos por un lado en SUMA1, y coger los denominadores y sumarlos por otro en SUMA2. 
#Despues imprime SUMA1 y SUMA2 en el fichero out_suma_fracciones_SNPs.xls

open (OUT, ">", 'out_suma_fracciones_SNPs.xls') or die "No puedo abrir out_suma_fracciones_SNPs.xls.\n";

print OUT "Suma1\tSuma2\n";

for (my $j=0;$j<=$#primera;$j++)
{

  if ( (index($primera[$j],"/") != -1) && (index($segunda[$j],"/") != -1)  )#Solo cojo aquellas lineas del fichero que contienen las dos fracciones
  {
    $posicion1 = index($primera[$j],"/");
    $posicion2 = index($segunda[$j],"/");
        
    $suma1 = substr($primera[$j],0,$posicion1+1) + substr($segunda[$j],0,$posicion2+1);
    $suma2 = substr($primera[$j],$posicion1+1,length($primera[$j])) + substr($segunda[$j],$posicion2+1,length($segunda[$j]));
    
    print OUT $suma1."\t".$suma2."\n";
  }
}                                                               

close (OUT);#Cierro out_suma_fracciones_SNPs.xls