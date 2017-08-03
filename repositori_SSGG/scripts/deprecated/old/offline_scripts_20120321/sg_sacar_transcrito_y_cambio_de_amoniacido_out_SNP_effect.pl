#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script extrae el cromosoma, la posicion, el transcrito y el cambio de aminoácido del output del script SNP_effect.pl\n";
        print "\nCOMO SE USA: sg_sacar_transcrito_y_cambio_de_amoniacido_out_SNP_effect.pl <fichero_SNP_effect>\n";
        print "ejemplo: sg_sacar_transcrito_y_cambio_de_amoniacido_out_SNP_effect.pl out_SNP_effect.txt\n\n";
	print "INPUT: fichero resultante de la anotación con SNP_effect.\n\n";
	print "OUTPUT: fichero tabulado en el que la primera columna es el transcrito y la segunda el cambio de aminoácido.\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}



# Abrimos los archivos de entrada
open(ENSEMBL,"<",$ARGV[0]);

my $substitution;

# Cargamos en memoria el fichero de la anotación del ensembl
for(my $j=1; my $ensembl=<ENSEMBL>;$j++)
{
	chomp($ensembl);
	push(@anotacion_ensembl,$ensembl);
}

for($n = 0; $n<=$#anotacion_ensembl;$n++)
{
		my @line_ensembl = split(/\t/,$anotacion_ensembl[$n]);
		my @sub = split(/\//,$line_ensembl[7]);
		my @pos = split(/:/,$line_ensembl[1]);
		$substitution = $sub[0].$line_ensembl[6].$sub[1];
		print $pos[0],"\t",$pos[1],"\t",$line_ensembl[3],"\t",$line_ensembl[4],"\t",$substitution,"\n";
}

exit;
