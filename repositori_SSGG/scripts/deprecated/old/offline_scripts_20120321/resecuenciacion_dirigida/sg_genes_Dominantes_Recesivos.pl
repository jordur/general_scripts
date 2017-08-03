#!/usr/bin/perl -w

#use List::MoreUtils qw(uniq);

sub usage {
	print "\nEste script crea un 2 archivos tabulados con la relación de genes Dominantes y Recesivos de  cada gen junto con su descriptcion.\n";
        print "\nCOMO SE USA: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl <fichero> col11 col2 col5 \n";
        print "ejemplo: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl Exoma50_MuestraA.xls 11 3 5\n\n";
	print "Está pensado para realizar el análisis de herencia recesiva o dominante con los ficheros de entrega de resecuenciación dirigida.\n\n";
	print "INPUT0: Fichero, ejemplo fichero final excel de resecuenciacion\n";
	print "INPUT1:El numero de columna donde se encuentra el ID del ENSEMBL en el INPUT, ejm ENSG00000005889\n";
	print "INPUT2:El numero de columna donde se encuentra la coordenada del SNPs o Indels\n";
	print "INPUT3:El numero de columna donde se encuentra la descripcion del gen\n";
	print "INPUT4:El numero de columna donde se ecuentra el Gene name\n";
	print "INPUT5:El numero de columna donde se encuentra la informacion homo/heterocigoto\n";
	print "OUTPUT: fichero tabulado en el que la primera columna es el código del gen y la segunda, el número de variaciones encontradas por gen y su descripcion. No hay que darle un fichero de salida.\n";
        exit(1);
	}	
# Si sólo ejecutamos el script, se imprime las instrucciones de uso

if(scalar(@ARGV) < 6){
    usage();
}

#DEFINIMOS LA COLUMNA DONDE ESTAN EL ID DEL ENSEMBL DE LOS GENES EJM ENSG00000005889

$input=$ARGV[0];
$ensg=$ARGV[1]-1;
$coord=$ARGV[2]-1;
$description=$ARGV[3]-1;
$geneName=$ARGV[4]-1;
$genotyper=$ARGV[5]-1;

%hash=();
%hash2=();
%hashDescription=();
%hashLines=();
%hashScalar=();
%hashName=();
%hashGenotyper=();

$string="";

open($REC_HEADER,'>',"Recessive_Candiate_Genes.txt"); # >1
open($REC_LINES,'>',"Recessive_Candiate_Genes_Lines.txt"); # >1

open($DOM_HEADER,'>',"Dominant_Candiate_Genes.txt"); # >=1
open($DOM_LINES,'>',"Dominant_Candiate_Genes_Lines.txt"); # >=1

%identify=();
open(FILE,'<',$input);
while($line=<FILE>)
{
        chomp($line);
        @explode=split("\t",$line);
	$string=$explode[$ensg]."_".$explode[$coord];

	$explode[$genotyper]=lc($explode[$genotyper]);
	$explode[$ensg]=lc($explode[$ensg]);

	if(!$identify{$string})
	{
		if($explode[$ensg])
	        {
			if($explode[$ensg]=~/^ensg/)
			{
				if(!$hash{$explode[$ensg]})
				{
					$hash{$explode[$ensg]}=1;
					$hashDescription{$explode[$ensg]}=$explode[$description];#DESCRIPTION OF GENE
					$hashName{$explode[$ensg]}=$explode[$geneName];

					if($explode[$genotyper]=~/^homo/)
					{
						$hashGenotyper{$explode[$ensg]}=1;
					}
				}
				else
				{
					$hash{$explode[$ensg]}=$hash{$explode[$ensg]}+1;
					if($explode[$genotyper]=~/^homo/)
					{
						$hashGenotyper{$explode[$ensg]}=1;
					}
				}
			}
			else
			{
				print "Gene no considerado en el analisis\t$explode[$ensg]\n";
			}
		}
		$identify{$string}=1;
	}
	if($explode[$ensg])
	{
		if($explode[$ensg]=~/^ensg/)
		{
			if(!$hash2{$explode[$ensg]})
			{
				$hash2{$explode[$ensg]}=1;
				$hashLines{$explode[$ensg]}[$hash2{$explode[$ensg]}]=$line;
				$hashScalar{$explode[$ensg]}=$hash2{$explode[$ensg]};
			}
			else
			{
				$hash2{$explode[$ensg]}=$hash2{$explode[$ensg]}+1;
				$hashLines{$explode[$ensg]}[$hash2{$explode[$ensg]}]=$line;
				$hashScalar{$explode[$ensg]}=$hash2{$explode[$ensg]};
			}
		}
	}
}
foreach $keys (keys %hash)
{

	$keys2=uc($keys);
	
	print $DOM_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";

	for($i=0;$i<($hashScalar{$keys}+1) ;$i++)
	{
		if($hashLines{$keys}[$i])
		{
			if($hashLines{$keys}[$i]=~/\W/)
			{
				print $DOM_LINES "$hashLines{$keys}[$i]\n";
			}
		}
	}
	
	if($hash{$keys}>1)
	{
	        for($i=0;$i<($hashScalar{$keys}+1) ;$i++)
      	 	{
	      	         if($hashLines{$keys}[$i])
              		 {
           	   	        if($hashLines{$keys}[$i]=~/\W/)
                      		{
                 	              	print $REC_LINES "$hashLines{$keys}[$i]\n";
                        	}
			}
		}
		print $REC_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";
	}
	else
	{
		if($hashGenotyper{$keys})
		{
			if($hashGenotyper{$keys}==1)
			{
				for($i=0;$i<($hashScalar{$keys}+1);$i++)
				{
					if($hashLines{$keys}[$i])
					{
						if($hashLines{$keys}[$i]=~/\W/)
						{
							print $REC_LINES "$hashLines{$keys}[$i]\n";
						}
					}

				}
				print $REC_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";
			}
		}
	}
}
