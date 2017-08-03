#!/usr/bin/perl -w

sub usage {
	print "\nEste script selecciona las SNVs comunes entre SNVs anotadas con el API del ensembl y SNVs anotadas con SIFT.\n";
        print "\nCOMO SE USA: sg_seleccionar_SNVs_comunes_anotacion_API_ensembl_y_SIFT.pl <fichero_ensembl> <fichero_sift>\n";
        print "ejemplo: sg_seleccionar_SNVs_comunes_anotacion_API_ensembl_y_SIFT.pl todos_chr_no_descritos_ensembl.txt SNVs_out_SIFTtxt\n\n";
	print "INPUT: Dos ficheros: el primero el fichero resultante de la anotación del ensembl. El segundo fichero es el resultante de la anotación de la página de SIFT.\n\n";
	print "OUTPUT: fichero tabulado SNVs_comunes_ensembl_SIFT_todos_chr.txt con las siguientes columnas:\n1: nombre del cromosoma\n2: coordenada en el cromosoma\n3: base en la referencia\n4: base detectada\n5: calidad del consenso\n6: calidad del SNP\n7: calidad máxima de mapeo\n8: coverage.
\n\n";
        exit(1);
	}	

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos el archivo de salida
open(OUT,">","SNVs_comunes_ensembl_SIFT_todos_chr.txt");

# Abrimos los archivos de entrada
open(ENSEMBL,"<",$ARGV[0]);
open(SIFT,"<",$ARGV[1]);

# Introducimos el fichero de la anotación con SIFT en un array para separar la primera columna
for(my $i=1; my $sift=<SIFT>;$i++)
{
	chomp($sift);
	my @todas_columnas = split(/\t/,$sift);
	if($todas_columnas[0] ne "Coordinates")
	{
		push(@snv_sift,$todas_columnas[0]);
	}
}

# Cerramos el archivo de SIFT
close(SIFT);

# Cargamos en memoria el fichero de la anotación del ensembl
for(my $j=1; my $ensembl=<ENSEMBL>;$j++)
{
	chomp($ensembl);
	push(@snv_ensembl,$ensembl);
}

# Cerramos el ficheros del Ensembl
close(ENSEMBL);

for(my $l=0; $l<=$#snv_sift; $l++)
{
	my @primera_columna = split(/,/,$snv_sift[$l]);
	my $chr = "chr".$primera_columna[0];
	for(my $n=0; $n<=$#snv_ensembl;$n++)
	{
		my @chr_coord_ensembl = split(/\t/,$snv_ensembl[$n]);
		if(($chr eq $chr_coord_ensembl[0]) && ($primera_columna[1] eq $chr_coord_ensembl[1]))
		{
			print OUT $chr_coord_ensembl[0],"\t",$chr_coord_ensembl[1],"\t",$chr_coord_ensembl[2],"\t",$chr_coord_ensembl[3],"\t",$chr_coord_ensembl[5],"\t",$chr_coord_ensembl[6],"\t",$chr_coord_ensembl[7],"\t",$chr_coord_ensembl[8],"\n";
		}
	}
}

# Cerramos los ficheros de 
close(ENSEMBL);

# Cerramos los ficheros de salida
close(OUT);

exit;

