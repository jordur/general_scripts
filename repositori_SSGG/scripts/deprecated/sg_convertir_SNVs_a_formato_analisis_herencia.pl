#!/usr/bin/perl -w

sub usage {
        my $usage =<<END;

PARA QUE SIRVE: parsea los archivos de SNVs a formato indel para poder compararlos 

COMO SE USA: sg_convertir_SNVs_a_formato_indel_analisis_herencia.pl <fichero_SNVs>

INPUT: archivo tabulado con el siguiente formato:

1       53580506        C       S       ENSG00000162383 ENST00000371491 A119P   -       NON_SYNONYMOUS_CODING             TOLERATED       0.38    SLC1A7  solute carrier family 1 (glutamate transporter), member 7 [Source:HGNC Symbol             ;Acc:10945]      101     53      73		0,405000001	OMIM    Comments


OUTPUT: archivo tabulado con el siguiente formato:

1       53580506   53580506        1      C       S     SNV  ENSG00000162383 ENST00000371491     NON_SYNONYMOUS_CODING   A119P     -      TOLERATED       0.38	SLC1A7  solute carrier family 1 (glutamate transporter), member 7 [Source:HGNC Symbol;Acc:10945]	N/A      N/A        73	101     35	0,405000001	OMIM	Comments


END

        print $usage;
}


# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
	usage();
}

# Abrimos los archivos de entrada
open(SNVS,"<",$ARGV[0]);

while($snv=<SNVS>)
{
	chomp($snv);
	@split_snv=split(/\t/,$snv);
	print $split_snv[0],"\t",$split_snv[1],"\t",$split_snv[1],"\t1\t",$split_snv[2],"\t",$split_snv[3],"\tSNV\t",$split_snv[4],"\t",$split_snv[5],"\t",$split_snv[8],"\t",$split_snv[6],"\t",$split_snv[7],"\t",$split_snv[9],"\t",$split_snv[10],"\t",$split_snv[11],"\t",$split_snv[12],"\tN/A\tN/A\t",$split_snv[15],"\t",$split_snv[13],"\t",$split_snv[14],"\t",$split_snv[16],"\t",$split_snv[17],"\t",$split_snv[18],"\n";
}

exit;
