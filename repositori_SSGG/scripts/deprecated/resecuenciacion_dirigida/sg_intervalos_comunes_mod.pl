#!/usr/bin/perl -w
use File::ReadBackwards;
#use strict;

my $no_coverage=0;
my $coordenada=0;
my $linea_input=0;
my $lineas_con_coverage=0;
my @cols;
my @columnas;
#my $ultima_linea_fichero  = File::ReadBackwards->new($ARGV[0])->readline;
my $coordenada_escrita;
my $cubiertos=$ARGV[1];
#chomp ($ultima_linea_fichero);


open(MATRIZ,"<",$ARGV[0]);
open(OUT,">","tempo_tmp");
while (my $laslineas=<MATRIZ>)
{
	chomp($laslineas);
	#$linea_input=$linea_input+1;
	@cols= split ("\t",$laslineas);
	#Calcula el numero de columnas distintas de 0, solamente aquellas posiciones cubiertas en todas las muestras seran tenidas en cuenta. Empieza en la tercera columna, la primera es el chr,la segunda la posicion. 
	if ($cubiertos eq "covered")
	{ 
		for (my $i=2; $i<=$#cols; $i++)
		{
			if ($cols[$i] != 0)
			{
				$no_coverage=$no_coverage+1;
			}
		}
		if ($no_coverage == $#cols-1) 
                {
                        print OUT $laslineas,"\n";
                }
                $no_coverage=0;
	}
	elsif ($cubiertos eq "not_covered")
        {
                for (my $i=2; $i<=$#cols; $i++)
                {
                        if ($cols[$i] == 0)
                        {
                                $no_coverage=$no_coverage+1;
                        }
                }
		if ($no_coverage == $#cols-1) 
                {
                        print OUT $laslineas,"\n";
                }
                $no_coverage=0;
        }
}
close (MATRIZ);
close (OUT);
my $ultima_linea_fichero  = File::ReadBackwards->new("tempo_tmp")->readline;
chomp ($ultima_linea_fichero);
open (NUEVO,"<","tempo_tmp");
while (my $posiciones = <NUEVO>)
{ 
	chomp ($posiciones);
	@columnas = split ("\t",$posiciones);
	#Este if es solo para imprimir la primera linea del fichero de intervalos.
		$lineas_con_coverage=$lineas_con_coverage+1;
		if  ($lineas_con_coverage == 1)
		{
			print $columnas[0],"\t",$columnas[1],"\t";
			$coordenada=$columnas[1];
		#	$coordenada_escrita=$columnas[1];
		}
		else
		{
			#5 posibilidades: 1) diferencia entre la coordenada anterior y la siguiente es mayor de 1 y no coincide con la ultima linea del fichero, se termina el intervalo anterior y se genera uno nuevo; 2) diferencia entre la coordenada anterior y la siguiente es igual a 1 y si coincide con la ultima linea del fichero, se termina el intervalo; 3) diferencia entre la coordenada anterior y la siguiente es mayor de 1 y si coincide con la ultima linea del fichero, se termina el intervalo y se genera el ultimo intervalo correspondiente a esa coordenada; 4) diferencia entre la coordenada anterior y la siguiente es igual a 1 y no coincide con la ultima linea del fichero, se cambia el valor de la variable $coordenada, que guarda el valor de la posicion anterior; 5) 
 			if ((($columnas[1] - $coordenada) > 1) && ($ultima_linea_fichero ne $posiciones))
			{
				print "$coordenada\n";
				print $columnas[0],"\t",$columnas[1],"\t";
				$coordenada=$columnas[1];
		#		$coordenada_escrita=$columnas[1];
			}
			elsif ((($columnas[1] - $coordenada) == 1) && ($ultima_linea_fichero eq $posiciones))
			{
				print $columnas[1],"\n";
			}
			elsif ((($columnas[1] - $coordenada) > 1) && ($ultima_linea_fichero eq $posiciones))
			{
				print $coordenada,"\n";
				print $columnas[0],"\t",$columnas[1],"\t",$columnas[1],"\n";
			}
			elsif ((($columnas[1] - $coordenada) == 1) && ($ultima_linea_fichero ne $posiciones))
			{
				$coordenada=$columnas[1];
			}
		}	
	#Contador de columnas con coverage igual a 0.
}
close(NUEVO);
system ('rm tempo_tmp');
sub usage 
{
        print "\nCOMO SE USA: \n\n";
        print "INPUT: Matriz del coverage por posicion y por muestra\n";
        print "OUTPUT: Intervalos\n\n";
        exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

sub bases_cubiertas
{
	my $decision = shift;
	if ($decision eq "cubiertas")
	{
		return "!= 0";
	}
	elsif ($decision eq "no_cubiertas")
	{
		return "== 0";
	}
}
