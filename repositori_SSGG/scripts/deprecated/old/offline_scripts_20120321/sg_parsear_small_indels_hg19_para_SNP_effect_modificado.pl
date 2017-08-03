#! /usr/bin/perl -w

sub usage {
	my $usage =<<END;

Script para parsear el listado de small indels obtenidos con Bioscope 1.2 tras convertir las coordenadas
de la versión hg18 a la hg19 con convertir_hg18_a_hg19.conf (/data/results/Solid0065/Exoma-HapMap-FA104-FA287/HapMap/bioscope230810_indels/small_indels_pairing_2/output/smallindel/anotacion_indels).

- CÓMO SE USA: sg_parsear_small_indels_hg19_para_SNP_effect_modificado.pl hg19_small_indels

- INPUT: listado de small indels tras convertir las coordenadas a hg19.

ejemplo: chr10   5203945 5203945 deletion        1       C/      5193943 GACCAA/GACAA    REF,4   4       HOMOZYGOUS      0.0039

- OUTPUT: Fichero tabulado con el formato adecuado para el script SNP_effect_modificado.pl con las siguientes columnas:
chr / start / end / reference / genotype / lengt / HEM / Score / Coverage / Type
Hay que dar un fichero de salida.

ejemplo: 10      5203945 5203945 C       -       1       HOMOZYGOUS      0.0039  4       deletion


END
	print $usage;
}

if(scalar(@ARGV) == 0){usage();}


# Abrimos el archivo de los small_indels
open(INDELS,'<',$ARGV[0]);

while (my $indel=<INDELS>)
{
	my $reference_allele = shift;
	my $genotype_allele = shift;
	chomp($indel);
	my @ex = split(/\t/,$indel);
	$ex[0] =~ tr/chr//d;
	$ex[5] =~ s/\///gi;
	if($ex[3] eq "deletion")
        {
                $reference_allele = $ex[5];
		$genotype_allele = "-";
	}
	elsif($ex[3] eq "insertion_site")
	{
		$reference_allele = "-";
		$genotype_allele = $ex[5];
	}
	print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$reference_allele,"\t",$genotype_allele,"\t",$ex[4],"\t",$ex[10],"\t",$ex[11],"\t",$ex[9],"\t",$ex[3],"\n";
}
