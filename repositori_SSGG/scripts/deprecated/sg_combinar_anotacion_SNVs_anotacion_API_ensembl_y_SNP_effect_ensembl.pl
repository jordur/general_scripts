#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script combina la información de las anotaciones del script SNP_effect.pl y de la API del Ensembl hecho con el script sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl. Este script sólo añade las columnas del coverage y la calidad del SNP. Lo utilizamos la primera vez que analizamos el exoma. En las versiones posteriores no se debería utilizar porque con la anotación de la API del ensembl ya debería estar solucionado.\n";
        print "\nCOMO SE USA: sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SNP_effect_ensembl <fichero_ensembl> <fichero_SNP_effect>\n";
        print "ejemplo: sg_seleccionar_SNVs_comunes_anotacion_API_ensembl_y_SNP_effect_ensembl.pl todos_chr_no_descritos_ensembl.txt SNVs_out_SNP_effect.txt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero resultante de la anotación del ensembl. El segundo fichero es el resultante de la anotación con SNP_effect.\n\n";
	print "OUTPUT: fichero tabulado SNVs_anotados_ensembl_combinando_API_y_script_SNP_effect.txt con las siguientes columnas:Chr\tCoordinate\tReference\tGenotype\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_ID\tQ_SNV\tCoverage\tConservation\n\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","mezcla");
print OUT "Chr\tCoordinate\tReference\tGenotype\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_ID\tQ_SNV\tCoverage\tConservation\n";


# Abrimos los archivos de entrada
open(ENSEMBL,"<",$ARGV[0]);
open(EFFECT,"<",$ARGV[1]);

# Introducimos el fichero de la anotación con SIFT en un array para separar la primera columna
for(my $i=1; my $effect=<EFFECT>;$i++)
{
	chomp($effect);
#	$effect =~ tr/:/\t/d;
	push(@anotacion_effect,$effect);
}

# Cerramos el archivo de SIFT
close(EFFECT);

my $substitution;
my @sub;

# Cargamos en memoria el fichero de la anotación del ensembl
for(my $j=1; my $ensembl=<ENSEMBL>;$j++)
{
	chomp($ensembl);
	push(@anotacion_ensembl,$ensembl);
}

# Cerramos el fichero del Ensembl
close(ENSEMBL);

for($k = 0; $k <=$#anotacion_effect; $k++)
{
	$anotacion_effect[$k]=~tr/:/\t/d;
	my @line_effect = split(/\t/,$anotacion_effect[$k]);
	for($n = 0; $n<=$#anotacion_ensembl;$n++)
	{
		my @line_ensembl = split(/\t/,$anotacion_ensembl[$n]);
		$line_ensembl[0]=~tr/chr//d;
		if(($line_effect[1] eq $line_ensembl[0]) && ($line_effect[2] eq $line_ensembl[1]))
		{
			if($line_effect[5] eq "SYNONYMOUS_CODING")
			{
				$substitution = $line_effect[8].$line_effect[7].$line_effect[8];
			}
			elsif($line_effect[5] eq "NON_SYNONYMOUS_CODING")
			{
				@sub = split(/\//,$line_effect[8]);
				$substitution = $sub[0].$line_effect[7].$sub[1];
			}
			else
			{
				$substitution = $line_effect[8];
			}
				print OUT $line_effect[1],"\t",$line_ensembl[1],"\t",$line_ensembl[2],"\t",$line_ensembl[3],"\t",$line_effect[4],"\t",$substitution,"\t",$line_effect[9],"\t",$line_effect[5],"\tN/A\tN/A\t",$line_effect[3],"\t",$line_ensembl[6],"\t",$line_ensembl[8],"\t",$line_ensembl[12],"\n";
#"\t",$line_ensembl[9],"\t",$line_ensembl[10],"\n";
		}
	}
}

system("sort -u mezcla > SNVs_anotados_ensembl_combinando_API_y_script_SNP_effect.txt");
system("rm mezcla");

# Cerramos los ficheros de salida
close(OUT);

exit;
