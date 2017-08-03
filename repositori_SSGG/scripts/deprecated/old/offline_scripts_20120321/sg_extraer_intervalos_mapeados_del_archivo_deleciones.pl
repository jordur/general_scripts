#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los intervalos de la secuencia referencia que han mapeado. Es el contrario al fichero de deleciones. Está pensado para convertir el fichero out_coverage en un fichero de intervalos de secuencia mapeada \n COMO SE USA:\n";
    print "INPUT: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el inicio, la segunda la coordenada el final. Las demás coordenadas no se tienen en cuenta. Debe estar ordenado de menor a mayor.\nejemplo: 343432      343445	12\n";
    print "       2) número del cromosoma al que pertenece. ejemplo: 1. No utilizar chr1.\n";
    print "       2) longitud del cromosoma. ejemplo: 247249719\n";
    print "\nEjemplo de uso:\nsg_extraer_intervalos_mapeados_del_archivo_deleciones.pl chr1_out_coverage.txt 1 247249719 > zonas_mapeadas_chr1\n\nHay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n\n";
    print " OUTPUT: Fichero tabulado en el que la primera columna es el cromosoma (formato chr1), la segunda es el inicio del intervalo mapeado y, la tercera, es el final del intervalo\nejemplo: chr1    247199329       247199353\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

my $chr = @ARGV[1];

my $longitud = @ARGV[2];

open(FILE,"<",@ARGV[0]);
my @fichero;
my @fila;
for (my $i=1; my $lines=<FILE>; $i++)
{
  chomp($lines);
  push (@fichero,$lines);
}
close (FILE);

my $l = 0;

for (my $j=0; $j<=$#fichero;$j++)
{
        @fila = split (/\t/,$fichero[$j]);
#	print "$fila[0]\t$fila[1]\n";
#	print "j----$j\n";
	if($j eq 0)
	{
		if($fila[0] eq 1)
		{
			my $inicial_1 = $fila[0] + 1;
			print "chr$chr\t$inicial_1\t";
		}
		elsif($fila[0] eq "Inicio")
		{
			@fila_siguiente = split (/\t/,$fichero[$j + 1]);
			if($fila_siguiente[0] eq 1)
               		{
				my $inicial_1 = $fila_siguiente[1] + 1;
                       		print "chr$chr\t$inicial_1\t";
				$j++;
               		}
			else
                	{
                        	my $final_0 = $fila_siguiente[0] - 1;
                        	my $inicial_0 = $fila_siguiente[1] + 1;
				print "chr",$chr,"\t1\t",$final_0,"\n";
				print "chr$chr\t$inicial_0\t";
				$j++;
			}
		}
		else
		{
			my $final_0 = $fila[0] - 1;
			my $inicial_0 = $fila[1] + 1;
			print "chr",$chr,"\t1\t",$final_0,"\n";
			print "chr$chr\t$inicial_0\t";
		}
	}
#	elsif($j eq 1)
#        {
#                if($fila[0] eq 1)
#                {
#                        my $inicial_1 = $fila[0] + 1;
#                        print "chr$chr\t$inicial_1\t";
#                }
#                else
#                {
#                        my $final_0 = $fila[0] - 1;
#                        my $inicial_0 = $fila[1] + 1;
#                        print "chr",$chr,"\t1\t",$final_0,"\n";
#                        print "chr$chr\t$inicial_0\t";
#                }
#
#        }
	elsif($j eq $#fichero)
	{
		if($fila[1] eq $longitud)
                {
			my $final = $fila[0] - 1;
			print "$final\n";
                }
                else
                {
			my $final = $fila[0] - 1;
			my $inicial = $fila[1] + 1;
			print "$final\n";
			print "chr$chr\t$inicial\t$longitud\n";
                }
	}
	else
	{
		my $final = $fila[0] - 1;
		my $inicial = $fila[1] + 1;
		print "$final\n";
		print "chr$chr\t$inicial\t";	
	}
}

exit;
