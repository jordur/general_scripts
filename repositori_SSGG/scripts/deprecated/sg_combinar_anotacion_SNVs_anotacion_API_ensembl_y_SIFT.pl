#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script combina la información de las anotaciones de SIFT y de la API del Ensembl.\n";
        print "\nCOMO SE USA: sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SIFT.pl <fichero_ensembl> <fichero_sift>\n";
        print "ejemplo: sg_seleccionar_SNVs_comunes_anotacion_API_ensembl_y_SIFT.pl todos_chr_no_descritos_ensembl.txt SNVs_out_SIFTtxt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero resultante de la anotación del ensembl. El segundo fichero es el resultante de la anotación de la página de SIFT.\n\n";
	print "OUTPUT: fichero tabulado SNVs_anotados_ensembl_SIFT_todos_chr.txt con las siguientes columnas:Chr\tCoordinate\tReference\tGenotype\tCodons\tTranscript_ID\tProtein_ID\tSubstitution\tRegion\tSNP_ID\tSNP_Type\tPrediction\tScore\tMedian_Info\t#Seqs_at_position\tGene_ID\tGene_Name\tGene_Description\tProtein_Family_ID\tProtein_Family_Description\tTranscript_Status\tProtein_Family_Size\tKa/Ks(Mouse)\tKa/Ks(Macaque)\tOMIM_Disease\tQ_SNV\tCoverage\tStrand\n";
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","anotacion");
print OUT "Chr\tCoordinate\tReference\tGenotype\tCodons\tTranscript_ID\tProtein_ID\tSubstitution\tRegion\tSNP_ID\tSNP_Type\tPrediction\tScore\tMedian_Info\t#Seqs_at_position\tGene_ID\tGene_Name\tGene_Description\tProtein_Family_ID\tProtein_Family_Description\tTranscript_Status\tProtein_Family_Size\tKa/Ks(Mouse)\tKa/Ks(Macaque)\tOMIM_Disease\tQ_SNV\tCoverage\tStrand\n";

# Codons\tSubstitution\tRegion\tSNP_ID\tSNV_Type\tPrediction\tScore\tMedian_info\tStrand\tTranscript_ID\tGene_ID\tGene_name\tGene_description\tProt_family_ID\tProt_family_description\tTranscript_status\tProt_family_size\t#Seqs_at_position\tConservation\tQ_consensus\tQ_SNV\tQ_max_map\tCoverage\tOMIM_Disease\n";

# Abrimos los archivos de entrada
open(ENSEMBL,"<",$ARGV[0]);
open(SIFT,"<",$ARGV[1]);

# Introducimos el fichero de la anotación con SIFT en un array para separar la primera columna
#for(my $i=1; my $sift=<SIFT>;$i++)
#{
#	chomp($sift);
#	push(@snv_sift,$sift);
#}

my @line_sift;
my $sub_split;
my %chr;

# Cargamos en memoria el fichero de la anotación del ensembl
for(my $j=1; my $ensembl=<ENSEMBL>;$j++)
{
	chomp($ensembl);
	push(@snv_ensembl,$ensembl);
}

# Cerramos el fichero del Ensembl
close(ENSEMBL);

while($snv_sift = <SIFT>)
{
	my $cont = 0;
	chomp($snv_sift);
#	$snv_sift=~ tr/,/\t/d;
	@line_sift = split(/\t/,$snv_sift);
	@sub_split = split(/,/,$line_sift[0]);
	$chr = "chr".$sub_split[0];
	for(my $n=0; $n<=$#snv_ensembl;$n++)
	{
		my @line_ensembl = split(/\t/,$snv_ensembl[$n]);
		if(($chr eq $line_ensembl[0]) && ($sub_split[1] eq $line_ensembl[1]))
		{
#			if($line_sift[16] eq $line_ensembl[9])
#			{
#				print "emparejados...",$chr,"=",$line_ensembl[0],".....",$line_sift[1],"=",$line_ensembl[1],"\n";
				print OUT $sub_split[0],"\t",$line_ensembl[1],"\t",$line_ensembl[2],"\t",$line_ensembl[3],"\t",$line_sift[1],"\t",$line_sift[2],"\t",$line_sift[3],"\t",$line_sift[4],"\t",$line_sift[5],"\t",$line_sift[6],"\t",$line_sift[7],"\t",$line_sift[8],"\t",$line_sift[9],"\t",$line_sift[10],"\t",$line_sift[11],"\t",$line_sift[12],"\t",$line_sift[13],"\t",$line_sift[14],"\t",$line_sift[15],"\t",$line_sift[16],"\t",$line_sift[17],"\t",$line_sift[18],"\t",$line_sift[19],"\t",$line_sift[20],"\t",$line_sift[21],"\t",$line_ensembl[6],"\t",$line_ensembl[8],"\t",$line_ensembl[12],"\n";
#			}
#			elsif ($line_sift[16] ne $line_ensembl[9])
#			{
#				if($line_sift[16] eq "")
#				{
#					print OUT $line_ensembl[0],"\t",$line_ensembl[1],"\t",$line_ensembl[2],"\t",$line_ensembl[3],"\t",$line_sift[4],"\t",$line_sift[7],"\t",$line_sift[8],"\t",$line_sift[9],"\t",$line_sift[10],"\t",$line_sift[11],"\t",$line_sift[12],"\t",$line_ensembl[11],"\t",$line_sift[5],"\t",$line_sift[15],"\t",$line_ensembl[9],"\t",$line_sift[13],"\t",$line_sift[17],"\t",$line_sift[18],"\t",$line_sift[19],"\t",$line_sift[20],"\t",$line_sift[21],"\t",$line_sift[14],"\t",$line_ensembl[12],"\t",$line_ensembl[5],"\t",$line_ensembl[6],"\t",$line_ensembl[7],"\t",$line_ensembl[8],"\t",$line_sift[24],"\n";
#					print "emparejados...",$chr,"=",$line_ensembl[0],".....",$line_sift[1],"=",$line_ensembl[1],".....",$line_sift[16],"=",$line_ensembl[9],".....nombre del gen a imprimir.....",$line_ensembl[9],"\n";
#				}
#				else
#				{
#					print "emparejados...",$chr,"=",$line_ensembl[0],".....",$line_sift[1],"=",$line_ensembl[1],".....",$line_sift[16],"=",$line_ensembl[9],".....nombre del gen a imprimir.....",$line_sift[16],"\n";
#				}
#			}
		$cont++;
		}if($cont>0){last;}
	}
}

system("sort -u anotacion > SNVs_anotados_ensembl_SIFT_todos_chr.txt");
system("rm anotacion");

# Cerramos los ficheros de salida
close(OUT);

exit;
