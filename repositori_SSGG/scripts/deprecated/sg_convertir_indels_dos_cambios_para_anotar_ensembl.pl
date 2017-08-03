#!/usr/bin/perl -w

sub usage {
        my $usage =<<END;

PARA QUE SIRVE: imprime los indels con dos cambios en líneas diferentes.

COMO SE USA: sg_convertir_indels_dos_cambios_para_anotar_ensembl.pl <fichero dos cambio>

INPUT: archivo tabulado con el siguiente formato:

chr	type		start		end		size	reference allele	sample allele
2       deletion        45645686        45645686        1       AAAAAAA 		GGGGGG/AAAAAA

OUTPUT: archivo tabulado para anotarlo con el script ensembl_variant_effect_predictor_v1-1.pl, con el siguiente formato:

2 	45645686	45645686	AAAAAAA/GGGGGG
2	45645686	45645686	AAAAAAA/AAAAAA

END

        print $usage;
}


# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos los archivos de entrada
open(INDELS,"<",$ARGV[0]);

while($indel=<INDELS>)
{
	chomp($indel);
	@split_indel=split(/\t/,$indel);
	$reference = "-";
	if($split_indel[5] =~ /^[ACGTN-]+$/)
	{
		$reference = $split_indel[5];
	}
	@split_sample=split(/\//,$split_indel[6]);
	if($split_indel[1] eq "insertion_site")
	{
		$start = $split_indel[2] + 1;
		print $split_indel[0],"\t",$start,"\t",$split_indel[3],"\t",$reference,"/",$split_sample[0],"\n";
		print $split_indel[0],"\t",$start,"\t",$split_indel[3],"\t",$reference,"/",$split_sample[1],"\n";
	}
	elsif($split_indel[1] eq "deletion")
	{
		print $split_indel[0],"\t",$split_indel[2],"\t",$split_indel[3],"\t",$reference,"/",$split_sample[0],"\n";
		print $split_indel[0],"\t",$split_indel[2],"\t",$split_indel[3],"\t",$reference,"/",$split_sample[1],"\n";
	}
	else
	{
		die("incorrect input");
	}
}

exit;
