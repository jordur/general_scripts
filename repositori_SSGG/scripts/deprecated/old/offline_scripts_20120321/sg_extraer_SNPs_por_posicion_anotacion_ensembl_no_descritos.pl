#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los SNPs por posición del fichero tabulado resultante de la anotación con la API del ensembl.\n\nCOMO SE USA:\nsg_extraer_SNPs_por_posicion_anotacion_ensembl_no_descritos.pl SNPs_no_descritos.txt\n\n";
    print " Input: 1) fichero tabulado de SNPs en el que la primera columna es el cromosoma, la segunda posición, la tercera la base referencia; la cuarta, la base encontrada; (...); la columna doce, la calidad del consenso; la trece, calidad del SNP; la catorce, la calidad máxima de mapeo; y la quince, el coverage. Debe estar ordenado de menor a mayor tanto en cromosoma como en posición.\nejemplo: chr1    17297200        C       Y       SYNONYMOUS_CODING       ENST00000375541 N/A     CROCC   ciliary rootlet coiled-coil, rootletin [Source:HGNC Symbol;Acc:21299]   ENST00000375541g.17297200Y>C    2.21000003814697        106     106     37      35\n\n";
    print " Output: SNPs_por_posicion_ensembl_no_descritos.txt. Es un fichero tabulado con las posiciones únicas de los SNPs no descritos en Ensembl. El resultado es un fichero tabulado donde la primera columna es el cromosoma, la segunda posición, la tercera la base referencia; la cuarta, la base encontrada; la quinta, la calidad del consenso; la sexta, calidad del SNP; la séptima, la calidad máxima de mapeo; y la octava, el coverage..\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(FILE,"<",@ARGV[0]);
my @fichero;
my @fila;

for (my $j=0;$lines=<FILE>;$j++)
{
	chomp($lines);
	push (@fichero,$lines);

}
close (FILE);

open(OUT,">","SNPs_por_posicion_ensembl_no_descritos.txt");

$k = 0;

for (my $i=0;$i<=$#fichero;$i++)
{
        @fila = split (/\t/,$fichero[$i]);

	if($fila[1] ne "position")
	{
		if($k ne $fila[1])
		{
			print OUT $fila[0],"\t",$fila[1],"\t",$fila[2],"\t",$fila[3],"\t",$fila[11],"\t",$fila[12],"\t",$fila[13],"\t",$fila[14],"\n";
			$k = $fila[1];
		}
		else 
		{
			$k = $fila[1];
		}
	}
}	

close(OUT);

exit;
