#! /usr/bin/perl -w

sub usage {
	my $usage =<<END;

Script para parsear el listado de small indels obtenidos con Bioscope 1.2 tras convertir las coordenadas
de la versión hg18 a la hg19 con convertir_hg18_a_hg19.conf (/data/results/Solid0065/Exoma-HapMap-FA104-FA287/HapMap/bioscope230810_indels/small_indels_pairing_2/output/smallindel/anotacion_indels).

- CÓMO SE USA: sg_parsear_small_indels_hg19_para_SNP_effect_modificado.pl hg19_small_indels

- INPUT: listado de small indels tras convertir las coordenadas a hg19.

ejemplo: chr10   5203945 5203945   C/   deletion 

- OUTPUT: Fichero tabulado con el formato adecuado para el script SNP_effect_modificado.pl con las siguientes columnas:
chr | start | end | reference/genotype | strand
Para crear las columnas nos basamos en el formato de entrada de la interfaz web del SNP_effect.pl(http://www.ensembl.org/Homo_sapiens/UserData/UploadVariations). Hay que dar un fichero de salida.

ejemplo: 10      5203945 5203945 C/-       1

END
	print $usage;
}

if(scalar(@ARGV) == 0){usage();}

my $reference_allele;
my $genotype_allele;
my @ex;

# Abrimos el archivo de los small_indels
open(INDELS,'<',$ARGV[0]);

while (my $indel=<INDELS>)
{
	chomp($indel);
	@ex = split(/\t/,$indel);
	$ex[0] =~ tr/chr//d;
	@split_ex3 = split(/\//,$ex[3]);
	if($split_ex3[0] eq "")
	{
		$split_ex3[0] = "-";
	}
	if($split_ex3[1] eq "")
        {
                $split_ex3[1] = "-";
        }
#	print "0...",$ex[0],"....5...",$ex[5],"....split3_0...",$split_ex3[0],"...split3_1...",$split_ex3[1],"\n";
#	print $ex[0],"\t",$ex[2],"\t",$ex[3],"\t",$split_ex3[0],"/",$split_ex3[1],"\t1\n";
	if($ex[4] eq "deletion")
        {
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$split_ex3[0],"/",$split_ex3[1],"\t1\n";
	}
	elsif($ex[4] eq "insertion_site")
	{
		print $ex[0],"\t",$ex[1]+1,"\t",$ex[2],"\t",$split_ex3[0],"/",$split_ex3[1],"\t1\n";
	}
#	print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$reference_allele,"\t",$genotype_allele,"\t1\n";
}

close(INDELS);


