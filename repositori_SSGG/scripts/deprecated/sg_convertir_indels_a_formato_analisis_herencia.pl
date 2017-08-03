#!/usr/bin/perl -w

sub usage {
        my $usage =<<END;

PARA QUE SIRVE: parsea los archivos de SNVs a formato indel para poder compararlos 

COMO SE USA: sg_convertir_SNVs_a_formato_indel_analisis_herencia.pl <fichero_SNVs>

INPUT: archivo tabulado con el siguiente formato:

1	10707901	10707901	1	-	G	insertion_site	ENSG00000130940	ENST00000377022	FRAMESHIFT_CODING	3771-3772	1151-1152	-	-	CASZ1	castor zinc finger 1 [Source:HGNC Symbol;Acc:26002]	G/GG/NO_CALL	REF,2,1		3	Comments	4,389999866

OUTPUT: archivo tabulado con el siguiente formato:

1	10707901	10707901	1	-	G	insertion_site	ENSG00000130940	ENST00000377022	FRAMESHIFT_CODING	-	-	N/A	N/A	CASZ1	Zinc finger protein castor homolog 1 (Castor-related protein)(Zinc finger protein 693) [Source:UniProtKB/Swiss-Prot;Acc:Q86V15]	G/GG/NO_CALL	REF,2,1		3	N/A	N/A	4,389999866	OMIM	Comments
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
	print $split_indel[0],"\t",$split_indel[1],"\t",$split_indel[2],"\t",$split_indel[3],"\t",$split_indel[4],"\t",$split_indel[5],"\t",$split_indel[6],"\t",$split_indel[7],"\t",$split_indel[8],"\t",$split_indel[9],"\t",$split_indel[12],"\t",$split_indel[13],"\tN/A\tN/A\t",$split_indel[14],"\t",$split_indel[15],"\t",$split_indel[17],"\t",$split_indel[18],"\t",$split_indel[19],"\tN/A\tN/A\t",$split_indel[21],"\t-\t",$split_indel[20],"\n";
}

exit;
