#!/usr/bin/perl -w

sub usage 
{
        print "PARA QUE SIRVE: Este programa es parte del pipeline para generar una matriz de posiciones cubierta/no cubiertas para un conjunto de muestras de manera que puedan identificarse la
s zonas que siempre se cubren entre ellas. El objetivo final es identificar las mutaciones que se localizan en las zonas que nunca se cubren con el kit de enriquecimiento que se este utilizando
\n";
        print "COMO SE USA: existe un wrapper que se llama sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango_wrapper.sh que incluye este y otro script \n";
        print "sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango_wrapper.sh fichero_pileup1 ... fichero_pileupN \n";
        print "La salida de este script es una matriz de ceros y unos para cada muestra y por cromosoma. Posteriormente, esta matriz pasa a formar parte de una matriz mayor en la que se incluye
n los resultados por muestra y por cromosoma asi como las posiciones\n";
        exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

my $contador=0;
my $str="basesNoCubiertas";
my $si=0;
my $sufijo = "rangos";

for(my $j = 0; $j <= $#ARGV; $j++)
{

	@miarray=();
        @inicio=();
        @final=();
        @split_primera_linea=();
        @ultima_linea=();

	#Abre fichero de salida	
	open (SALIDA,">",$ARGV[$j].$str);
	#Abre fichero de input, matrices de 1s y 0s atendiendo a si una base esta presente o no en el fichero pileup
	open(MATRIZ,"<",$ARGV[$j]);
		$cromosoma = $ARGV[$j];
		$cromosoma =~s /.*chr//;
		$cromosoma =~s /\.tab//;
		chomp ($cromosoma);
        	while (my $posiciones=<MATRIZ>)
        	{
			chomp ($posiciones);
			@split_linea = split ("\t",$posiciones);
			for (my $i =2; $i<= $#split_linea; $i++)
			{
				if ($split_linea[$i] eq "1")
				{
					$contador = $contador+1;
				}
			}
			$porcentaje = ($contador*100/($#split_linea-1));
			#print $porcentaje,"<-------PORCENTAJEEEE\n";
		#if ($contador != ($#split_linea-1))
		if ($porcentaje < 80)
			{
				print SALIDA $split_linea[0],"\t",$split_linea[1],"\n";
			}
			$contador=0;

		}
	close (MATRIZ);
	close (SALIDA);
	open (ENTRA,"<",$ARGV[$j].$str);
	open (SALGO,">",$ARGV[$j].$sufijo);
	while (my $lineas = <ENTRA>)
	{
		chomp ($lineas);
		push (@miarray,$lineas);
	}
	if (scalar (@miarray) < 1)
	{
		last;
	}
	
	@split_primera_linea= split ("\t",$miarray[0]);
	push (@inicio,$split_primera_linea[1]);
	for (my $nu = 0; $nu < $#miarray; $nu++)
	{
		@posiciones_antes = split ("\t",$miarray[$nu]);
		@posiciones_dp = split ("\t",$miarray[$nu+1]);
		if (($posiciones_dp[1]-$posiciones_antes[1])>1)
        	{
			push(@final,$posiciones_antes[1]);#Guardo las coordenadas del final de la deleccion
            		push(@inicio,$posiciones_dp[1]);#Guardo las coordenadas del inicio de la deleccion
        	}
	}		
	@ultima_linea= split ("\t",$miarray[-1]);
	push(@final,$ultima_linea[1]);
	for (my $pa=0;$pa<=$#inicio;$pa++)
	{
		#print  "ESTO ES FICHERO",$ARGV[$j].$sufijo ,"chr".$cromosoma."_+_".$inicio[$pa]."_".$final[$pa]."\t",$inicio[$pa],"\t",$final[$pa],"\t+\n";
		print SALGO  "chr".$cromosoma."_".$inicio[$pa]."_".$final[$pa]."_+\t",$inicio[$pa],"\t",$final[$pa],"\t+\n";
	}
	close (ENTRA);
	close (SALGO);
}
