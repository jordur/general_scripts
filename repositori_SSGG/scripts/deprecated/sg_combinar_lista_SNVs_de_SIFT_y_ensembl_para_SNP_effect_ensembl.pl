#!/usr/bin/perl -w

sub usage {
	print "\nEste script crea un archivo en formato el formato para anotarlo con el script SNP_effect.pl desarrollado por el ensembl.\n";
        print "\nCOMO SE USA: sg_combinar_lista_SNVs_de_SIFT_y_ensembl_para_SNP_effect_ensembl.pl <fichero_ensembl> <fichero_sift>\n";
        print "ejemplo: sg_combinar_lista_SNVs_de_SIFT_y_ensembl_para_SNP_effect_ensembl.pl SNVs_ensembl.txt SNVs_sift.txt\n\n";
	print "INPUT: Dos ficheros. El primero es el resultante de realizar la anotación con el script sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl. El segundo es el fichero resultante de la anotación con SIFT.\n\n";
	print "OUTPUT: fichero con el formato 2 36871	36871	T/G 	1	\nEl listado se recoge en el fichero SNVs_para_SNP_effect.txt con el que podemos ejecutar directamente SNP_effect.\n\n";
        exit(1);

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
}
if(scalar(@ARGV) == 0){
    usage();
}

# Abrimos el archivo de salida
open(OUTPUT,">","combinacion");

# Abrimos el archivo del SIFT
open(SIFT,"<",$ARGV[1]);

for(my $i =1;my $sift = <SIFT>;$i++)
{
	chomp $sift;
	$sift =~ tr/,/\t/d;
	push(@lineas_sift,$sift);
}
close(SIFT);

# Abrimos el archivo del ensembl
open(ENSEMBL,"<",$ARGV[0]);

while (my $ensembl=<ENSEMBL>)
{
        chomp($ensembl);
	my @ex = split(/\t/,$ensembl);
	$ex[0] =~ tr/chr//d;
	for(my $j=0; $j<=$#lineas_sift; $j++)
	{
		my @ex_sift = split(/\t/,$lineas_sift[$j]);
		if(($ex[0] eq $ex_sift[0]) && ($ex[1] eq $ex_sift[1]))
		{
		print OUTPUT "$ex[0]\t$ex[1]\t$ex[1]\t$ex[2]/$ex[3]\t$ex[11]\n";
		}
	} 
}

close(ENSEMBL);

# Cerramos los ficheros de salida
close(OUTPUT);

system("sort -u combinacion > SNVs_para_SNP_effect.txt");
system("rm combinacion");

exit;
