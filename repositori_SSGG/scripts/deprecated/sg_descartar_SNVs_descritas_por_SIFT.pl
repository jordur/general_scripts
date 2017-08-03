#!/usr/bin/perl -w

sub usage 
{
	print "\nEste script combina la información de las anotaciones de SIFT y de la API del Ensembl.\n";
        print "\nCOMO SE USA: sg_descartar_SNVs_descritas_por_SIFT.pl <fichero_ensembl> <fichero_sift>\n";
        print "ejemplo: sg_descartar_SNVs_descritas_por_SIFT.pl todos_chr_no_descritos_ensembl.txt SNVs_out_SIFTtxt\n\n";
	print "INPUT: Dos ficheros: el primero es la combinación de la acotación del SNP_effect.pl y el script sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl para combinar la anotación del API del ensembl.Ejemplo:\n10      100152870       G       C       ENST00000408492 -       -       DOWNSTREAM      N/A     N/A     ENSG00000221419 hsa-mir-1287 hsa-mir-1287 [Source:miRBase;Acc:MI0006349]     -1      18      3       0.634999990463257\nEl segundo fichero es el resultante de la anotación de la página de SIFT.Ejemplo:\n10,103991381,1,G/T      -       ENST00000370002 ENSP00000359019 NA      EXON CDS        rs2281983:A     Unknown Not scored  NA       NA      NA      ENSG00000107859 PITX3   Pituitary homeobox 3 (Paired-like homeodomain transcription factor 3)(Homeobox protein PITX3) [Source:UniProtKB/Swiss-Prot;Acc:O75364]       ENSFM00260000050475     PITUITARY HOMEOBOX PAIRED  HOMEODOMAIN TRANSCRIPTION FACTOR HOMEOBOX KNOWN   7       0.007   0.015    CATARACT, POSTERIOR POLAR, 4; PAIRED-LIKE HOMEODOMAIN TRANSCRIPTION FACTOR 3; ANTERIOR SEGMENT MESENCHYMAL DYSGENESIS       :       :\n\n
Importante: El fichero de la anotación del SIFT no debe tener líneas repetidas!!!!\n";
	print "OUTPUT: Dos ficheros: \n- SNVs_anotados_ensembl_SIFT_todos_chr_no_descritos.txt: SNVs no anotadas como conocidas por el SIFT, con el mismo formato de entrada que la combinación de la anotación del Ensembl.\n- SNVs_anotados_ensembl_SIFT_todos_chr_descritos.txt: SNVs anotadas como conocidas por el SIFT, con el mismo formato de entrada que la combinación de la anotación del Ensembl.\n"; 
        exit(1);
}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","anotacion");

open(DESCRITOS,">","anotacion_descritos");

#open(SIFT_MODIF,">","sift");

# Abrimos los archivos de entrada
open(ENSEMBL,"<",$ARGV[0]);
open(SIFT,"<",$ARGV[1]);
#while(my $linea_sift=<SIFT>)
#{
#        chomp($linea_sift);
#	$linea_sift=~ tr/,/\t/d;
#	my @split_linea_sift = split(/\t/,$linea_sift);
#	$split_linea_sift[0] =~s/chr//gi;
#	print SIFT_MODIF $split_linea_sift[0],"\t",$split_linea_sift[1],"\t",$split_linea_sift[9],"\n";
#}

#close(SIFT_MODIF);
#system("sort -u sift > sift_unicos");

#open(SIFT_MODIF2,"<","sift_unicos");

# Cargamos en memoria el fichero de la anotación del ensembl
for(my $j=1; my $sift=<SIFT>;$j++)
{
	chomp($sift);
	$sift=~ tr/,/\t/d;
	push(@snv_sift,$sift);
}

# Cerramos el fichero del Ensembl
close(SIFT);

my $chr;
my @line_sift;
my @split_ensembl;
my $cont = 0;

while($snv_ensembl = <ENSEMBL>)
#while($snv_sift = <SIFT_MODIF2>)
#for(my $t=0; $t<=$#nuevo_sift;$t++)
{
	chomp($snv_ensembl);
#	my $chr = "chr".$line_sift[0];
	for(my $n=0; $n<=$#snv_sift;$n++)
	{
		@line_sift = split(/\t/,$snv_sift[$n]);
		$chr = $line_sift[0];
		@line_ensembl = split(/\t/,$snv_ensembl);
		if(($chr eq $line_ensembl[0]) && ($line_sift[1] eq $line_ensembl[1]))
		{
#			if(($line_sift[9] eq "novel") || ($line_sift[9] eq "NA"))
#			{	
#				print OUT $line_ensembl[0],"\t",$line_ensembl[1],"\t",$line_ensembl[2],"\t",$line_ensembl[3],"\t",$line_ensembl[4],"\t",$line_ensembl[5],"\t",$line_sift[9],"\t",$line_ensembl[7],"\t",$line_ensembl[8],"\t",$line_ensembl[9],"\t",$line_ensembl[10],"\t",$line_ensembl[11],"\t",$line_ensembl[12],"\t",$line_ensembl[13],"\n";
			$cont++;
#			}
#			else
#			{
			print DESCRITOS $line_ensembl[0],"\t",$line_ensembl[1],"\t",$line_ensembl[2],"\t",$line_ensembl[3],"\t",$line_ensembl[4],"\t",$line_ensembl[5],"\t",$line_sift[9],"\t",$line_ensembl[7],"\t",$line_ensembl[8],"\t",$line_ensembl[9],"\t",$line_ensembl[10],"\t",$line_ensembl[11],"\t",$line_ensembl[12],"\t",$line_ensembl[13],"\t",$line_ensembl[14],"\t",$line_ensembl[15],"\t-\n";
#             		}
		}
	}
	if($cont == 0)
	{
#		print "cont no descritos.....$cont\n";
		print OUT $snv_ensembl,"\n";
	}
	elsif($cont > 0)
	{
		$cont = 0;
	}
}

system("sort -u anotacion > SNVs_anotados_ensembl_SIFT_todos_chr_no_descritos.txt");
system("rm anotacion");
system("sort -u anotacion_descritos > SNVs_anotados_ensembl_SIFT_todos_chr_descritos.txt");
system("rm anotacion_descritos");

# Cerramos los ficheros de salida
close(OUT);
close(DESCRITOS);

exit;
