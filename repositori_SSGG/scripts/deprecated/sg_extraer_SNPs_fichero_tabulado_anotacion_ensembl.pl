#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los SNPs de un fichero tabulado que se encuentran en las regiones indicadas por un archivo donde indicamos los intervalos de secuencias que estamos analizando \n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) fichero tabulado de SNPs en el que la primera columna es el cromosoma, la segunda posición, la tercera la base referencia; la cuarta, la base encontrada; y a partir de la quinta toda la información correspondiente a la anotación de la API del ensembl\n";
    print " Output: fichero tabulado con los SNPs que se encuentran en las regiones indicadas en el fichero de input 1. Hay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FILE1,"<",@ARGV[0]);
open(FILE2,"<",@ARGV[1]);
my @fichero1;
my @fichero2;
my @fila1;
my @fila2;
for (my $k=1; my $lines=<FILE1>; $k++)
{
  chomp($lines);
  push (@fichero1,$lines);
}

close (FILE1);


for (my $n=1; my $lines2=<FILE2>; $n++)
{
  chomp($lines2);
  push (@fichero2,$lines2);
}

close (FILE2);

my $l = 0;

for (my $i=0; $i<=$#fichero1;$i++)
{
        @fila1 = split (/\t/,$fichero1[$i]);
#	print "fila1...j=$j....i=$i....l=$l....",$fila1[0],"\t",$fila1[1],"\t",$fila1[2],"\n";
	for (my $j=0; $j<=$#fichero2;)
	{
		@fila2 = split (/\t/,$fichero2[$l]);
#		$j = $l;
#		print "fila2..j=$j....i=$i....l=$l...",$fila2[0],"\t",$fila2[1],"\t",$fila2[2],"\t",$fila2[3],"\n";
		if($fila1[1] > $fila2[1])
                {
		$j++;
		$l++;
                }
		elsif(($fila1[1] <= $fila2[1]) && ($fila1[2] >= $fila2[1]))
		{
			print $fila2[0],"\t",$fila2[1],"\t",$fila2[2],"\t",$fila2[3],"\t",$fila2[4],"\t",$fila2[5],"\t",$fila2[6],"\t",$fila2[7],"\t",$fila2[8],"\t",$fila2[9],"\t",$fila2[10],"\t",$fila2[11],"\t",$fila2[12],"\t",$fila2[13],"\t",$fila2[14],"\n";
			$j++;
			$l++;
		}
		elsif($fila2[1] > $fila1[2])
		{	
			last;
		}
	}
}

exit;
