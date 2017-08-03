#!/usr/bin/perl
sub usage {
    print "\nPARA QUE SIRVE: Chequea si existen lecturas que mapean entre dos coordenadas\n";
    print "Input: 1) fichero .csfasta con las lecturas mapeadas , 2) rango inicio , 3) rango final\n";
    print "Output: Número de lecturas en ese rango \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(MACSFASTA,"<",@ARGV[0]);
my $inicio = @ARGV[1];
my $final = @ARGV[2];
my $contador=0;
my $contador_lineas =0;
while (my $lineas = <MACSFASTA>)
{
	chomp ($lineas);
	$primer_caracter = substr ($lineas, 0, 1);
	@split_posicion = split (",",$lineas);
	$contador_lineas = $contador_lineas+1;	
	#print "linea numero", $contador_lineas,"\n";
	if ($primer_caracter eq ">")
	{
		#print "holaaa\n";
		for (my $i =1; $i<=$#split_posicion; $i++)
		{
			#print $#split_posicion,"\n"; # $split_posicion[$i],"\n";
			$split_posicion[$i] =~s/.*_//;
			$split_posicion[$i] =~s/\.[0-9]//;
			$split_posicion[$i] =~s/-//;
			#print $split_posicion[$i],"\n";
			#$calcular_posicion = &posicion($split_posicion[$i]);
			if (($split_posicion[$i] >= $inicio) && ($split_posicion[$i] <= $final))
			{
				$contador=$contador+1;
				#print "contador parcial",$contador,"\n";
#				push (@lecturas_en_esas_coordenadas, $lineas);
			}
		}	
	}
}

print $contador,"\n";
#my @split_posicion2;
#sub posicion ($pos)
#{
#	@split_posicion2 = split ("_",$pos);
#	$split_posicion2[1] =~s/\.*$//;
#	$split_posicion2[1] =~s/-//;
#	return $split_posicion2[1];

#}
#my $contador=0;

