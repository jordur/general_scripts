#!/usr/bin/perl -w
use File::ReadBackwards;

sub usage 
{
	print "PARA QUE SIRVE: Este programa es parte del pipeline para generar una matriz de posiciones cubierta/no cubiertas para un conjunto de muestras de manera que puedan identificarse las zonas que siempre se cubren entre ellas. El objetivo final es identificar las mutaciones que se localizan en las zonas que nunca se cubren con el kit de enriquecimiento que se este utilizando\n";
	print "COMO SE USA: existe un wrapper que se llama sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango_wrapper.sh que incluye este y otro script \n";
	print "sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango_wrapper.sh fichero_pileup1 ... fichero_pileupN \n";
	print "La salida de este script es una matriz de ceros y unos para cada muestra y por cromosoma. Posteriormente, esta matriz pasa a formar parte de una matriz mayor en la que se incluyen los resultados por muestra y por cromosoma asi como las posiciones\n";
        exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}

open(POSICIONES,$ARGV[$#ARGV]);
	
	while (my $posiciones=<POSICIONES>)
	{
		chomp($posiciones);
		push (@posiciones_referencia,$posiciones);
	}

close(POSICIONES);


$str = "comp";
for($j = 0; $j <= $#ARGV-1; $j++)
{
	
	open (SALIDA,">",$ARGV[$j].$str);
	#print SALIDA join ("\n",@posiciones_referencia),"\n";
	open (FICHERO,"<",$ARGV[$j]) or die ("ack - $!");
	# print "----ABRO FICHERO----\n";
	my $ultima_linea_fichero  = File::ReadBackwards->new($ARGV[$j])->readline;
	if (!$ultima_linea_fichero)
	{
	#	print "FICHERO NO EXISTE\n";
		for ($p=0; $p <= $#posiciones_referencia ; $p++)
		{
			print SALIDA "0\n";#$posiciones_referencia[$p],"--0\n";
		}	
	}
	else
	{
	chomp ($ultima_linea_fichero);
	my $val=0;
	while ($posiciones_muestra = <FICHERO>)
	{
		chomp ($posiciones_muestra);
		for ($mo=$val; $mo <= $#posiciones_referencia ; $mo++)
		{
			chomp ($posiciones_referencia[$mo]);
	
#		print "LO QUE ENTRA-REF-MUESTRA--> ",$posiciones_referencia[$mo],"\t",$posiciones_muestra,"\n";
			if($posiciones_referencia[$mo] eq $posiciones_muestra)
			{
				if (($mo < $#posiciones_referencia) && ($posiciones_referencia[$mo] eq $ultima_linea_fichero))
				{
					$resto=$#posiciones_referencia-$mo;
					for ($quedan =0; $quedan <=$resto; $quedan++)
					{
	#					print "ULTIMAS-->",$posiciones_referencia[$mo+$quedan],"--0\n";
						print SALIDA "0\n";
					}
					last;
					
				}
				else
				{
	#				print "ES IGUAL-->",$posiciones_referencia[$mo],"--1\n";
					print SALIDA "1\n";
					$val=$mo+1;
					last;
				}
			}
			elsif ($posiciones_muestra > $posiciones_referencia[$mo])
			{
	#			print "MUESTRA MAYOR IMPRIMO LA BASE----> ",$posiciones_referencia[$mo],"--0\n";
				print SALIDA "0\n";
			}
		}
	}
	}
	close(FICHERO);
	close(SALIDA);
}	
