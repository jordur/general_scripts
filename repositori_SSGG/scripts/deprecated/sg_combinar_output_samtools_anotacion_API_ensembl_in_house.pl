#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script combina la información de las anotaciones del script hecho in house de la anotacion del ensembl (sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl) con la información original de las SNVs detectadas y convertidas al hg19. Este script sólo añade las cuatro columnas.\n";
        print "\nCOMO SE USA: sg_combinar_output_samtools_anotacion_API_ensembl_in_house.pl <fichero_samtools_hg19> <fichero_anotacion_API_in_house>\n";
        print "ejemplo: sg_combinar_output_samtools_anotacion_SNP_effect.pl todas_SNVs_hag19_formato_samtools.txt snp_no_descritos_todos_chr.txt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero del samtools con las coordenadas en la misma versión que la anotación con la API. El segundo fichero es el resultante de la anotación con sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl.\n\n";
	print "OUTPUT: fichero tabulado SNVs_descritas_ensembl_combinando_samtools_y_anotacion_API_ensembl_in_house.txt con las siguientes columnas:chr\tcoordinate\treference\tgenotype\tSNP_ID\tQ_consensus\tQ_SNV\tQ_max_map\tcoverage\tgene_name\tgene_descrition\tstrand\tconservation\n\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","SNVs_descritas_ensembl_combinando_samtools_y_anotacion_API_ensembl_in_house.txt");
# Sin dar formato a la substitución del aminoácido en los "NON_SYNONYMOUS_CODING" y "SYNONYMOUS_CODING"
# print OUT "Chr\tCoordinate\tReference\tGenotype\tGene_ID\tTranscript_ID\tConsequence\tPosition in cDNA\tPosition in protein\tAmino acid change\tSNP_ID\tQ_SNV\tCoverage\n";
# Dando formato a la substitución del aminoácido en los "NON_SYNONYMOUS_CODING" y "SYNONYMOUS_CODING"
print OUT "chr\tcoordinate\treference\tgenotype\tSNP_ID\tQ_consensus\tQ_SNV\tQ_max_map\tcoverage\tgene_name\tgene_descrition\tstrand\tconservation\n";

# Abrimos los archivos de entrada
open(SAMTOOLS,"<",$ARGV[0]);
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

# Cargamos en memoria el fichero del samtools
for(my $j=1; my $samtools=<SAMTOOLS>;$j++)
{
	chomp($samtools);
	push(@resultado_samtools,$samtools);
}

# Cerramos el fichero del Ensembl
close(SAMTOOLS);

for($k = 0; $k <=$#anotacion_effect; $k++)
{
	@line_effect = split(/\t/,$anotacion_effect[$k]);
	for($n = 0; $n<=$#resultado_samtools;$n++)
	{
		@line_samtools = split(/\t/,$resultado_samtools[$n]);
		if(($line_effect[0] eq $line_samtools[0]) && ($line_effect[1] eq $line_samtools[1]))
		{
			print OUT $line_effect[0],"\t",$line_effect[1],"\t",$line_effect[2],"\t",$line_samtools[3],"\t",$line_effect[4],"\t",$line_samtools[4],"\t",$line_samtools[5],"\t",$line_samtools[6],"\t",$line_samtools[7],"\t",$line_effect[9],"\t",$line_effect[10],"\t",$line_effect[11],"\t",$line_effect[12],"\n";
		}
	}
}

#system("sort -u mezcla > SNVs_descritas_ensembl_combinando_samtools_y_script_SNP_effect.txt");
#system("rm mezcla");

# Cerramos los ficheros de salida
close(OUT);

exit;
