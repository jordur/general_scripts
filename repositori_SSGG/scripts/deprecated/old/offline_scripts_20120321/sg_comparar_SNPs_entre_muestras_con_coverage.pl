#!/usr/bin/perl -w

# Este script tiene como salida un archivo en que cada línea es una posición del cromosoma donde se ha localizado un
# SNP en alguna muestra. El archivo resultante es matriz_snps.txt, en la que la segunda columna son todas las 
# posiciones del cromosoma donde se ha encontrado un SNP. Las restantes columnas son los SNPs encontrados en cada
# muestra. Si no se ha encontrado un SNP en la posición, se imprime un cero. 
# El script crea un fichero temporal /tmp que se borra al final de la ejecución.


sub usage {
        print "\nCOMO SE USA: sg_comparar_SNPs_entre_muestras_con_coverage.pl <archivo_1> <archivo_2> ... <archivo_N>\n\n";
        print "ejemplo: g_comparar_SNPs_entre_muestras_con_coverage.pl snp1.txt snp2.txt snp3.txt\n\n";
	print "INPUT: ficheros con los SNPs detectados en las muestras que queramos comparar. Deben ser ficheros tabulados cuya primera posición sea el cromosoma, la segunda la posición, la tercera, la base en la secuencia referencia; la cuarta, el genotipo; la quinta, la calidad del consenso; la sexta, calidad del SNP; la séptima la calidad máxima de mapeo; y la octava, el coverage.\n\n";
	print "OUTPUT: matriz_snps.txt\n\n";
	print "La salida del programa es un fichero en el que encontramos todas las posiciones de todos los SNPs encontrados en todos los ficheros. La primera columna son todas las coordenadas. La segunda, la base referencia; A partir de la tercera columna están las muestras. De cada muestra se imprimen dos columnas, la base encontrada y el coverage en la posición del SNPs. Si en la muestra encontramos el SNP, imprime la base del genotipo, sinó, imprime cero.\n\n";
        exit(1);
}

# Si no introducimos el fichero de entrada, se imprime el modo de uso
if(scalar(@ARGV) == 0){
    usage();
}

####################################################################################
#
# Creamos el fichero que contiene las posiciones de todos los SNPs y lo cargamos en memoria
#
####################################################################################

# Creamos el directorio temporal
system("mkdir tmp");

# Extraemos la segunda columna de todos los archivos en el directorio temporal
for($i = 0; $i <= $#ARGV; $i++)
{
	system("awk '{print \$2\"\t\"\$3}' $ARGV[$i] > tmp/posiciones_$ARGV[$i]");
}

# Unimos todos los ficheros, los ordenamos y nos quedamos con las filas únicas
system("cat tmp/posiciones_* | sort -u -k1n > tmp/posiciones.txt");

# Abrimos el archivo de las posiciones y lo cargamos en memoria
open(POSICIONES,"tmp/posiciones.txt");
	
while (my $posiciones=<POSICIONES>)
{
	chomp($posiciones);
	push (@todas_posiciones,$posiciones);
}

close(POSICIONES);


######################################################################################
#
# Creamos los ficheros de una columna donde imprimimos las bases del genotipo para cada posición.
#
######################################################################################

# Arbrimos archivo a archivo
for($j = 0; $j <= $#ARGV; $j++)
{
	# Abrimos cada fichero mediante mediantela subrutina abrir_fichero
	@posiciones_snp = abrir_fichero($ARGV[$j]);

	# Abrimos los ficheros donde vamos a imprimir los resultados
	open(MUESTRA,">","tmp/snps_para_matriz_$ARGV[$j]");
#	print MUESTRA $ARGV[$j],"\n";
	
	$m = 0;	
	$n = 0; # la variable $n nos va a marcar los SNPs que no coinciden con las posiciones.txt
	# Comparamos cada posición de SNP frente a la de posiciones.txt
	for ($k= 0; $k <= $#todas_posiciones;$k++)
	{
		@posiciones_referencia = split(/\t/,$todas_posiciones[$k]);
#		$m = 0;
		for ($l= $m; $l <= $#posiciones_snp;)
		{
			# Introducimos cada línea/SNP en un array
			@posiciones_muestra = split(/\t/,$posiciones_snp[$l]);
			
				# Si la posición de la referencia coincide con nuestro SNP, imprimimos la base
				if($posiciones_referencia[0] eq $posiciones_muestra[1])
				{
					print MUESTRA $posiciones_muestra[3],"\t",$posiciones_muestra[7],"\n";
#					print $posiciones_referencia[0],"-----------",$posiciones_muestra[1],"--------base:$posiciones_muestra[3]\n";
#					print $posiciones_muestra[3],"\n";
					last;
				}
				# Si la posición de la referencia no coincide con nuestro SNP, imprime el valor 0
				elsif($posiciones_referencia[0] < $posiciones_muestra[1])
				{
					print MUESTRA "0\t0\n";
#					print $posiciones_referencia[0],"-----------",$posiciones_muestra[1],"----------base: 0\n";
#					print "0\n";
#					$l = $m;	# Si no hay un SNP, igualamos los valores de estas variables para que empiece con el siguiente valor y no salte valores
					$n++;
					last;
				}
				# Esta condición no sirve para nada. Lee los valores pero no imprime nada
				elsif($posiciones_referencia[0] > $posiciones_muestra[1])
				{	
#					 print $posiciones_referencia[0],"-----------",$posiciones_muestra[1],"--------referencia mayor\n";		
					$l++;
				}
				$m++;
#			}
		}

	}
	# Completamos el fichero con tantos ceros como posiciones de referencia quedan
	$numero_snp = $#posiciones_snp + 1;
	$finales = $k - $n - $numero_snp;
	if($k > $n + $numero_snp)
	{
		for($o = 1; $o <= $finales; $o++)
		{
#			print "0\n";
			print MUESTRA "0\t0\n";
#print "--------------------------------------------------------------------------------------------------------\n";
		}
	}

	# Cerramos el fichero de salida
	close(MUESTRA);
}

# Pegamos todos los ficheros en uno llamado matriz_snps.txt
system("paste tmp/posiciones.txt tmp/snps_para_matriz_* > matriz_snps.txt");
system("rm -rf tmp");


# Salimos del programa
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
