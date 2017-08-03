#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script completa la anotación realizada con la API del ensembl y SNP_effect.pl en los que no está descrito el external name que utiliza ensembl. Es el script complementario al sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SNP_effect_ensembl.pl si no tenemos la información de la cadena. Para utilizarlo, debemos tener una lista de los genes y transcritos en los que hemos encontrado las variaciones.\n";
        print "\nCOMO SE USA: sg_agregar_informacion_genes_usando_biomart.pl <fichero_biomart> <fichero_SNP_anot>\n";
        print "ejemplo: sg_agregar_informacion_genes_usando_biomart.pl genes_biomart.txt SNVs_out_SNP_anot.txt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero tabulado resultante de la anotación del biomart con las columnas ID_gen_Ensembl | ID_transcript_Ensembl | Gene_external_name | Strand | Gen_Description \n Ejemplo: ENSG00000001561 ENST00000371401 ENPP4   1       ectonucleotide pyrophosphatase/phosphodiesterase 4 (putative function) [Source:HGNC Symbol;Acc:3359] \nEl segundo fichero es el resultante de la anotación con el script sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SNP_effect_ensembl.pl.\n\n";
	print "OUTPUT: fichero tabulado SNVs_anotados_biomart_SIFT_todos_chr.txt con las siguientes columnas:Chr\tCoordinate\tReference\tGenotype\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_ID\tGene_Name\tGene_Description\tStrand\tQ_SNV\tCoverage\tConservation\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
#open(OUT,">","remezcla");
open(OUT,">","SNVs_anotados_biomart_combinando_API_y_script_SNP_anot_con_genes.txt");
print OUT "Chr\tCoordinate\tReference\tGenotype\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_ID\tGene_Name\tGene_Description\tStrand\tQ_SNV\tCoverage\tConservation\n";


# Abrimos los archivos de entrada
open(BIOMART,"<",$ARGV[0]);
open(ANOT,"<",$ARGV[1]);

# Introducimos el fichero de la anotación con SIFT en un array para separar la primera columna
for(my $i=1; my $anot=<ANOT>;$i++)
{
	chomp($anot);
#	$anot =~ tr/:/\t/d;
	push(@anotacion_anot,$anot);
}

# Cerramos el archivo de SIFT
close(ANOT);


# Cargamos en memoria el fichero de la anotación del biomart
for(my $j=1; my $biomart=<BIOMART>;$j++)
{
	chomp($biomart);
	push(@anotacion_biomart,$biomart);
}

# Cerramos el fichero del Ensembl
close(BIOMART);

for($k = 0; $k <=$#anotacion_anot; $k++)
{
	my @line_anot = split(/\t/,$anotacion_anot[$k]);
	for($n = 0; $n<=$#anotacion_biomart;$n++)
	{
		my @line_biomart = split(/\t/,$anotacion_biomart[$n]);
		if(($line_anot[4] eq $line_biomart[1]) && ($line_anot[10] eq $line_biomart[0]))
		{
				print OUT $line_anot[0],"\t",$line_anot[1],"\t",$line_anot[2],"\t",$line_anot[3],"\t",$line_anot[4],"\t",$line_anot[5],"\t",$line_anot[6],"\t",$line_anot[7],"\t",$line_anot[8],"\t",$line_anot[9],"\t",$line_anot[10],"\t",$line_biomart[2],"\t",$line_biomart[4],"\t",$line_biomart[3],"\t",$line_anot[11],"\t",$line_anot[12],"\t",$line_anot[13],"\n";
		}
	}
}

#system("sort -u remezcla > SNVs_anotados_biomart_combinando_API_y_script_SNP_anot_con_genes.txt");
#system("rm remezcla");

# Cerramos los ficheros de salida
close(OUT);

exit;
