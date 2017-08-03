#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Combina la anotación de las SNVs filtradas por SIFT y realizada con ensembl_variant_effect_predictor_v1-1.pl, combinada con samtools y biomart con la anotación de SIFT para evitar la anotación repetida. Para ello, comparamos los valores de posición y transcrito. Sólo sirve para combinar las variaciones no sinónimas y no descritas en el dbSNP provenientes de SIFT\n COMO SE USA:\nsg_unir_anotacion_samtools_SNP_effect_biomart_SIFT.pl <fichero_SIFT_no_sinonimas_unicas> <fichero_ensembl_variation_samtools_biomart>\n
ejemplo: sg_unir_anotacion_samtools_SNP_effect_biomart_SIFT.pl unicos_SIFT_no_sinonimos.txt SNVs_anotados_ensembl_SIFT_todos_chr_no_descritos.txt\n";
    print "Input: 1) fichero tabulado con las SNVs no sinónimas anotadas con SIFT y combinada con la anotación de la API del Ensembl (generada con el script sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SIFT.pl)\nejemplo:\n
10,102743792,1,A/T      TGG-aGG ENST00000342071 ENSP00000339844 W169R   EXON CDS        novel   Nonsynonymous   TOLERATED0.19     4.32    2       ENSG00000055950 MRPL43  39S ribosomal protein L43, mitochondrial Precursor (L43mt)(MRP-L43)(Mitochondrial ribosomal protein bMRP36a) [Source:UniProtKB/Swiss-Prot;Acc:Q8N983]     ENSFM00250000003611     39S RIBOSOMAL L43 MITOCHONDRIAL PRECURSOR L43MT MRP L43 MITOCHONDRIAL RIBOSOMAL BMRP36A   KNOWN   11      0.211   0.311\n";
    print "       2) fichero tabulado con las SNVs anotadas con ensembl_variant_effect_predictor_v1-1.pl, combinada con samtools y biomart\nejemplo\n
10      102056809       G       K       ENST00000318222 I371I   -       SYNONYMOUS_CODING       N/A     N/A     ENSG00000193      PKD2L1  polycystic kidney disease 2-like 1 [Source:HGNC Symbol;Acc:9011]        -1      60      110     -
\n";
    print " Output: fichero tabulado con las SNVs anotadas y sin repeticiones. El fichero de salida se llama SNVs_no_descritas_anotadas_final.txt\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


system("sort -nk2 @ARGV[1] > ensembl_ordenado");

open(SIFT_ORIGINAL,"<","@ARGV[0]");
open(ENSEMBL,"<","ensembl_ordenado");
open(INTER,">","inter");
open(OUT,">","combinacion");
print OUT "Chr\tCoordinate\tReference allele\tSample allele\tGene_ID\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_Name\tGene_Description\tQ_SNV\tCoverage\tGERP Conservation\tOMIM_disease\n"; 

my @sift;
my @ensembl;

my @split_sift;
my @split_ensembl;

for (my $k=1; my $line_sift=<SIFT_ORIGINAL>; $k++)
{
	chomp($line_sift);
	@split_sift=split(/\t/,$line_sift);
	@split_sift_primera=split(/,/,$split_sift[0]);
#Sólo imprimimos de SIFT lo que nos interesa para comparar e imprimir: chr, coord, ref/genotype, trancript_ID, SNP_type, Prediction, Score, OMIM.
	@sift_elegido=(@split_sift_primera[0,1],@split_sift[2,7,8,9,21]);
	print INTER join("\t",@sift_elegido),"\n";
}

close (SIFT_ORIGINAL);

system("sort -nk2 inter > sift_ordenado");

open(SIFT,"<","sift_ordenado");
for (my $hu=1; my $line_sift_def=<SIFT>; $hu++)
{
	chomp($line_sift_def);
	push(@sift,$line_sift_def);
	
}
close(SIFT);

for (my $n=1; my $line_ensembl=<ENSEMBL>; $n++)
{
	chomp($line_ensembl);
	push (@ensembl,$line_ensembl);
}

close (ENSEMBL);

my $snp_type;
my $l=0;
my $h=0;
my $b = 0;
my $o = 0;
my $p = 0;

for(my $i=0; $i<=$#sift;$i++)
{
	@split_sift=split (/\t/,$sift[$i]);
	for(my $j=0; $j<=$#ensembl;)
	{
		$j=$l;
		@split_ensembl=split(/\t/,$ensembl[$l]);
		if($split_sift[1] > $split_ensembl[1])
		{
			print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[10],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
			$l++;
			$j++;
                }
		elsif($split_sift[1] < $split_ensembl[1])
		{
			last;
		}
		elsif ($split_ensembl[1] == $split_sift[1])
		{
			if($split_ensembl[4] eq $split_sift[2])
			{
				if((($split_sift[3] eq "Nonsynonymous") && ($split_ensembl[7] eq "NON_SYNONYMOUS_CODING")) || (($split_sift[3] eq "Nonsynonymous") && ($split_ensembl[7] eq "STOP_GAINED")) || (($split_sift[3] eq "Nonsynonymous") && ($split_ensembl[7] eq "STOP_LOST")) || (($split_sift[3] eq "Nonsynonymous") && ($split_ensembl[7] eq "NMD_TRANSCRIPT")))
                                {
					print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[10],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\tSIFT-",$split_ensembl[7],"\t",$split_sift[4],"\t",$split_sift[5],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\t",$split_sift[6],"\n";
					$l++;
                                        $j++;
                                }
                                else
                                {
					print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[10],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
					$l++;
                        		$j++;
				}
			}
			else
			{	
			print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[10],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
			$l++;
			$j++;
			}
		}
	}	
	
}

if($l<$#ensembl)
{
	for(my $y=$l; $y<=$#ensembl;$y++)
	{
		@split_ensembl=split(/\t/,$ensembl[$y]);
		print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[10],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
	}
}

my $com = " '/^\$/d'";

close(OUT);

system("sed $com combinacion | sort -u > SNVs_no_descritas_anotadas.txt");
system("rm ensembl_ordenado sift_ordenado combinacion inter");

exit;
