#!/usr/bin/perl -w


sub usage {
        print "\nCOMO SE USA: sg_comparar_SNPs_entre_muestras.pl <archivo_1> <archivo_2> ... <archivo_N>\n\n";
	print "INPUT: ficheros con SNPs o el fichero pileup de las muestras que se quieren comparar. Deben ser ficheros tabulados cuya primera posición sea el cromosoma, la segunda la posición, la tercera la base en la secuencia referencia y la cuarta, el genotipo. A partir de la quinta columna, el contenido no interviene en el script.\n\n";
	print "OUTPUT: matriz_snps.txt\n\n";
	print "La salida del programa es un fichero en el que encontramos todas las posiciones de todos los SNPs encontrados en todos los ficheros. La primera columna son todas las coordenadas. A partir de la segunda columna están las muestras. Si en la muestra encontramos el SNP, imprime la base del genotipo, sinó, imprime cero.\n\n";
        exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}


system("mkdir tmp");

# Extraemos la segunda columna de todos los archivos en el directorio temporal
for($i = 0; $i <= $#ARGV; $i++)
{
	system("awk '{print \$2}' $ARGV[$i] > tmp/posiciones_$ARGV[$i]");
}

# Unimos todos los ficheros, los ordenamos y nos quedamos con las filas únicas
system("cat tmp/posiciones_* | sort -u -n > tmp/posiciones.txt");

# Abrimos el archivo de las posiciones y lo cargamos en memoria
open(POSICIONES,"tmp/posiciones.txt");
	
while (my $posiciones=<POSICIONES>)
{
	chomp($posiciones);
	push (@todas_posiciones,$posiciones);
}

close(POSICIONES);



for($j = 0; $j <= $#ARGV; $j++)
{
	@posiciones_snp = abrir_fichero($ARGV[$j]);

	open(MUESTRA,">","tmp/snps_para_matriz_$ARGV[$j]");
	
	$m = 0;	
	$n = 0; 
	for ($k= 0; $k <= $#todas_posiciones;$k++)
	{
		@posiciones_referencia = split(/\t/,$todas_posiciones[$k]);
		for ($l= $m; $l <= $#posiciones_snp;)
		{
			@posiciones_muestra = split(/\t/,$posiciones_snp[$l]);
			
				# Si la posición de la referencia coincide con nuestro SNP, imprimimos la base
				if($posiciones_referencia[0] eq $posiciones_muestra[1])
				{
    					print MUESTRA $posiciones_muestra[3],"\n";
					last;
				}
				# Si la posición de la referencia no coincide con nuestro SNP, imprime el valor 0
				elsif($posiciones_referencia[0] < $posiciones_muestra[1])
				{
					print MUESTRA "0\n";
					$n++;
					last;
				}
				elsif($posiciones_referencia[0] > $posiciones_muestra[1])
				{	
					$l++;
				}
				$m++;
		}

	}
	# Si el numero de posiciones es mayor que el numero de posicion en el fichero pileup 
	$numero_snp = $#posiciones_snp + 1;
	$finales = $k - $n - $numero_snp;
	if($k > $n + $numero_snp)
	{
		for($o = 1; $o <= $finales; $o++)
		{
			print MUESTRA "0\n";
		}
	}

	close(MUESTRA);
}

system("paste tmp/posiciones.txt tmp/snps_para_matriz_* > matriz_snps.txt");
system("rm -rf tmp");


exit;

##################################################################################################
#
#     SUBRUTINAS
#
##################################################################################################

sub abrir_fichero
{
        my @snp_leidos;
        open(SNP,"<",@_);
        while (my $snp=<SNP>)
        {
        chomp($snp);
        push (@snp_leidos,$snp);
        }
        close(SNP);
        return @snp_leidos;
}

##########################################
