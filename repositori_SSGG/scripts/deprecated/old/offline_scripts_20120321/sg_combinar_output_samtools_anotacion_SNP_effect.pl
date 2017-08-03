#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script combina la información de las anotaciones del script SNP_effect.pl y la información original de las SNVs detectadas (output de samtools). Este script sólo añade las columnas del coverage y la calidad del SNP.\n";
        print "\nCOMO SE USA: sg_combinar_output_samtools_anotacion_SNP_effect.pl <fichero_samtools> <fichero_SNP_effect>\n";
        print "ejemplo: sg_combinar_output_samtools_anotacion_SNP_effect.pl chr10SNVs_hg19 SNVs_out_SNP_effect.txt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero del samtools. El segundo fichero es el resultante de la anotación con SNP_effect.pl.\n\n";
	print "OUTPUT: fichero tabulado SNVs_descritas_ensembl_combinando_samtools_y_script_SNP_effect.txt con las siguientes columnas:Chr\tCoordinate\tReference\tGenotype\tGene_ID\tTranscript_ID\tConsequence\tPosition in cDNA\tSubstitution\tSNP_ID\tQ_consensus\tQ_SNV\tQ_max_map\tCoverage\n\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","mezcla");
# Sin dar formato a la substitución del aminoácido en los "NON_SYNONYMOUS_CODING" y "SYNONYMOUS_CODING"
# print OUT "Chr\tCoordinate\tReference\tGenotype\tGene_ID\tTranscript_ID\tConsequence\tPosition in cDNA\tPosition in protein\tAmino acid change\tSNP_ID\tQ_SNV\tCoverage\n";
# Dando formato a la substitución del aminoácido en los "NON_SYNONYMOUS_CODING" y "SYNONYMOUS_CODING"
print OUT "Chr\tCoordinate\tReference\tGenotype\tGene_ID\tTranscript_ID\tConsequence\tPosition_cDNA\tSubstitution\tSNP_ID\tQ_consensus\tQ_SNV\tQ_max_map\tCoverage\n";

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
	$anotacion_effect[$k]=~tr/:/\t/d;
	@line_effect = split(/\t/,$anotacion_effect[$k]);
	for($n = 0; $n<=$#resultado_samtools;$n++)
	{
		@line_samtools = split(/\t/,$resultado_samtools[$n]);
		$line_samtools[0]=~s/chr//gi;
		if(($line_effect[1] eq $line_samtools[0]) && ($line_effect[2] eq $line_samtools[1]))
		{

# Si no se quiere introducir el cambio aminoacídico para "NON_SYNONYMOUS_CODING", "STOP_GAINED", "STOP_LOST", "NMD_TRANSCRIPT", "SPLICE_SITE" y "SYNONYMOUS_CODING" eliminar 

# desde aquí...
			if($line_effect[5] eq "SYNONYMOUS_CODING")
			{
				$substitution = $line_effect[8].$line_effect[7].$line_effect[8];
			}
			elsif(($line_effect[5] eq "NON_SYNONYMOUS_CODING") || ($line_effect[5] eq "STOP_GAINED") || ($line_effect[5] eq "STOP_LOST"))
			{
				@sub = split(/\//,$line_effect[8]);
				$substitution = $sub[0].$line_effect[7].$sub[1];
			}
			elsif(($line_effect[5] eq "NMD_TRANSCRIPT") || ($line_effect[5] eq "SPLICE_SITE"))
			{
				if($line_effect[8] ne "-")
				{
					if(length($line_effect[8]) == 1)
					{
						$substitution = $line_effect[8].$line_effect[7].$line_effect[8];
					}
					else
					{
						@sub = split(/\//,$line_effect[8]);
						$substitution = $sub[0].$line_effect[7].$sub[1];
					}
				}
				else
				{
					$substitution = $line_effect[8];
				}
			}
			else
			{
				$substitution = $line_effect[8];
			}
			print OUT $line_samtools[0],"\t",$line_samtools[1],"\t",$line_samtools[2],"\t",$line_samtools[3],"\t",$line_effect[3],"\t",$line_effect[4],"\t",$line_effect[5],"\t",$line_effect[6],"\t",$substitution,"\t",$line_effect[9],"\t",$line_samtools[4],"\t",$line_samtools[5],"\t",$line_samtools[6],"\t",$line_samtools[7],"\n";

#  ... hasta aquí



# Sin dar formato a la substitución del aminoácido en los "NON_SYNONYMOUS_CODING" y "SYNONYMOUS_CODING" 
#
#			print OUT $line_effect[1],"\t",$line_samtools[1],"\t",$line_samtools[2],"\t",$line_samtools[3],"\t",$line_effect[3],"\t",$line_effect[4],"\t",$line_effect[5],"\t",$line_effect[6],"\t",$line_effect[7],"\t",$line_effect[8],"\t",$line_effect[9],"\t",$line_effect[10],"\t",$line_samtools[5],"\t",$line_samtools[7],"\n";
		}
	}
}

system("sort -u mezcla > SNVs_descritas_ensembl_combinando_samtools_y_script_SNP_effect.txt");
system("rm mezcla");

# Cerramos los ficheros de salida
close(OUT);

exit;
