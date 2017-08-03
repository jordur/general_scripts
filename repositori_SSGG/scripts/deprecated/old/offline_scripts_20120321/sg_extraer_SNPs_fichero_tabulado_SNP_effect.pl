#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los SNPs anotados con ensembl_variant_effect_predictor_v1-1.pl de un fichero tabulado que se encuentran en las regiones indicadas por un archivo donde indicamos los intervalos de secuencias que estamos analizando \n\nCOMO SE USA: sg_extraer_SNPs_fichero_tabulado_SNP_effect.pl <fichero_intervalos> <fichero_SNP_effect>\n\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) output del script de anotación del Ensembl ensembl_variant_effect_predictor_v1-1.pl. Es un fichero tabulado SNPs en el que la primera columna es el nombre de la variación, la segunda chr:posición, el gen, etc. Permite los dos formatos de entrada.Ejemplos:\n1_1158631_A/G   1:1158631       ENSG00000078808 ENST00000263741 SYNONYMOUS_CODING       863     190     D       -\n1_1158631_A/G   1	1158631       ENSG00000078808 ENST00000263741 SYNONYMOUS_CODING       863     190     D       -\n";
    print " Output: fichero tabulado con los SNPs que se encuentran en las regiones indicadas en el fichero de input 1. Hay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(REGION,"<",@ARGV[0]);
open(SNV,"<",@ARGV[1]);

# Definimos variables
my @intervalo;
my @snv;
my @split_intervalo;
my @split_snv;


for (my $n=1; my $lines2=<SNV>; $n++)
{
  chomp($lines2);
  push (@snv,$lines2);
}

close (SNV);


while($intervalo=<REGION>)
{
	chomp($intervalo);
        @split_intervalo = split (/\t/,$intervalo);
	for (my $j=0; $j<=$#snv;$j++)
	{
#		Sustituimos ":" por tabulación por si nos sirven directamente del output de ensembl_variant_effect_predictor_v1-1.pl o viene separado por cromosomas.
		$snv[$j]=~s/:/\t/g;
		@split_snv = split (/\t/,$snv[$j]);
#		@sub_split_snv = split(/:/, $split_snv[1]);
		if(($split_intervalo[1] <= $split_snv[2]) && ($split_intervalo[2] >= $split_snv[2]))
                {
			print $snv[$j],"\n";
		}
	}
}

close(REGION);

exit;
