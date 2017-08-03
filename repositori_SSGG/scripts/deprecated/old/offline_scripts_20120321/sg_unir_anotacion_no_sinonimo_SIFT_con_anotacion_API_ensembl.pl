#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Combina la anotación realizada con la API de Ensembl y la anotación de SIFT para evitar la anotación repetida. Para ello, comparamos los valores de posición y transcrito. Sólo sirve para combinar las variaciones no sinónimas y no descritas en el dbSNP provenientes de SIFT\n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las SNVs no sinónimas anotadas con SIFT y combinada con la anotación de la API del Ensembl (generada con el script sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SIFT.pl)\n\n";
    print "       2) fichero tabulado con las SNVs anotadas con la API del Ensembl (generado con el script sg_combinar_anotacion_SNVs_anotacion_API_ensembl_y_SNP_effect_ensembl.pl o con la segunda versión de la anotación con la API del Ensembl)\n";
    print " Output: fichero tabulado con las SNVs anotadas y sin repeticiones. El fichero de salida se llama SNVs_no_descritas_anotadas.txt\n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

system("sort -nk2 @ARGV[0] > sift_ordenado");
system("sort -nk2 @ARGV[1] > ensembl_ordenado");

open(SIFT,"<","sift_ordenado");
open(ENSEMBL,"<","ensembl_ordenado");
open(OUT,">","combinacion");
print OUT "Chr\tCoordinate\tReference Genotype\tTranscript_ID\tSubstitution\tSNP_ID\tSNP_Type\tPrediction\tScore\tGene_ID\tGene_Name\tGene_Description\tQ_SNV\tCoverage\tConservation\tOMIM_disease\n"; 

my @sift;
my @ensembl;

my @split_sift;
my @split_ensembl;

for (my $k=1; my $line_sift=<SIFT>; $k++)
{
	chomp($line_sift);
#	$line_sift =~ s/,/\t/gi;
	push (@sift,$line_sift);
}

close (SIFT);


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
	@split_sift = split (/\t/,$sift[$i]);
#	print $split_sift[1],"\n";
#	@split_sift_siguiente = split (/\t/,$sift[$p]);
	for(my $j=0; $j<=$#ensembl;)
	{
		$j=$l;
		@split_ensembl = split (/\t/,$ensembl[$l]);
#		$h = $l + 1;
#		@split_ensembl_siguiente = split (/\t/,$ensembl[$h]);
		if($split_sift[1] > $split_ensembl[1])
		{
#			print OUT $ensembl[$l],"\n";
			print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[10],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
			$l++;
			$j++;
                }
		elsif($split_sift[1] < $split_ensembl[1])
		{
			last;
		}
		elsif ($split_ensembl[1] == $split_sift[1])
		{
			if($split_ensembl[4] eq $split_sift[5])
			{
#				print "entra...sift...",$split_sift[1],"...ensembl...",$split_ensembl[1],"\n";
				my @split_sift2 = split (/\t/,$sift[$i]);
				if((($split_sift2[10] eq "Nonsynonymous") && ($split_ensembl[7] eq "NON_SYNONYMOUS_CODING")) || (($split_sift2[10] eq "Nonsynonymous") && ($split_ensembl[7] eq "STOP_GAINED")) || (($split_sift2[10] eq "Nonsynonymous") && ($split_ensembl[7] eq "STOP_LOST")) || (($split_sift2[10] eq "Nonsynonymous") && ($split_ensembl[7] eq "NMD_TRANSCRIPT")))
                                {
#					print "entra...",$split_sift2[10],"....",$split_ensembl[7],"\n";
                                        $snp_type = $split_sift2[10]."-".$split_ensembl[7];
					print OUT $split_sift2[0],"\t",$split_sift2[1],"\t",$split_sift2[2],"\t",$split_sift2[3],"\t",$split_sift2[5],"\t",$split_sift2[7],"\t",$split_sift2[9],"\t",$snp_type,"\t",$split_sift2[11],"\t",$split_sift2[12],"\t",$split_sift2[15],"\t",$split_sift2[16],"\t",$split_sift2[17],"\t",$split_sift2[25],"\t",$split_sift2[26],"\t",$split_sift2[27],"\t",$split_sift2[24],"\n";
					$l++;
                                        $j++;
                                }
                                else
                                {
#					print "entra2...",$split_sift2[10],"....",$split_ensembl[7],"\n";
					print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[10],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
					$l++;
                        		$j++;
				}
#                                        $snp_type = $split_sift2[8];
#					print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[10],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
 #                               }
#				$b++;
#				$l++;
#				last;
			}
			else
			{	
			print OUT $split_ensembl[0],"\t",$split_ensembl[1],"\t",$split_ensembl[2],"\t",$split_ensembl[3],"\t",$split_ensembl[4],"\t",$split_ensembl[5],"\t",$split_ensembl[6],"\t",$split_ensembl[7],"\t",$split_ensembl[8],"\t",$split_ensembl[9],"\t",$split_ensembl[10],"\t",$split_ensembl[11],"\t",$split_ensembl[12],"\t",$split_ensembl[14],"\t",$split_ensembl[15],"\t",$split_ensembl[16],"\n";
#			print OUT $ensembl[$l],"\n";
			$l++;
			$j++;
			}
		}
	}	
}

close (OUT);

my $com = " '/^\$/d'";

system("sed $com combinacion | sort -u > SNVs_no_descritas_anotadas.txt");
system("rm ensembl_ordenado sift_ordenado combinacion");

exit;
