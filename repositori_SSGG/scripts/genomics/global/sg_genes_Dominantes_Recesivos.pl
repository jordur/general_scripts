#!/usr/bin/perl -w

#use List::MoreUtils qw(uniq);
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use Exporter;

sub usage {
	print "\nEste script crea un 2 archivos tabulados con la relacion de genes Dominantes y Recesivos de cada gen junto con su descriptcion.\n";
	print "\nCOMO SE USA: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl <fichero> col11 col2 col5 \n";
	print "ejemplo: sg_extraer_genes_unicos_y_numero_SNV_por_gen.pl Exoma50_MuestraA.xls 11 3 5\n\n";
	print "Esta pensado para realizar el analisis de herencia recesiva o dominante con los ficheros de entrega de resecuenciacion dirigida.\n\n";
	print "INPUT 1: Fichero, ejemplo fichero final excel de resecuenciacion\n";
	print "INPUT 2: El numero de columna donde se encuentra el nombre del gen (HGNC o ENSEMBL) en el INPUT, ejm TNC1,ENSG00000005889\n";
	print "INPUT 3: El numero de columna donde se encuentra la coordenada del SNPs o Indels\n";
	print "INPUT 4: El numero de columna donde se encuentra la descripcion del gen\n";
	print "INPUT 5: El numero de columna donde se encuentra la informacion homo/heterocigoto\n";
	print "INPUT 6: Separador entre columnas del fichero de anotacion\n";
	print "INPUT 7: Prefijo para la salida\n";
	print "OUTPUT: fichero tabulado en el que la primera columna es el codigo del gen y la segunda, el numero de variaciones encontradas por gen y su descripcion. No hay que darle un fichero de salida.\n";
	exit(1);
	}	
# Si solo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) < 6){
    usage();
}

#DEFINIMOS LA COLUMNA DONDE ESTAN EL ID DEL ENSEMBL DE LOS GENES EJM ENSG00000005889

my $input=$ARGV[0];
my $ensg=$ARGV[1]-1;
my $coord=$ARGV[2]-1;
my $description=$ARGV[3]-1;
my $geneName=$ARGV[1]-1;
my $genotyper=$ARGV[4]-1;
my $sep=$ARGV[5];
my $prefix=$ARGV[6];

my %hash=();
my %hash2=();
my %hashDescription=();
my %hashLines=();
my %hashScalar=();
my %hashName=();
my %hashGenotyper=();

my $string="";

open(my $REC_HEADER,'>',$prefix."Recessive_Candidate_Genes.txt"); # >1
open(my $REC_LINES,'>',$prefix."Recessive_Candidate_Genes_Lines.txt"); # >1

open(my $DOM_HEADER,'>',$prefix."Dominant_Candidate_Genes.txt"); # >=1
open(my $DOM_LINES,'>',$prefix."Dominant_Candidate_Genes_Lines.txt"); # >=1

my %identify=();
open(FILE,'<',$input);
my $line_number = 0;
while(my $line=<FILE>) {
	$line_number ++;
	if ($line_number > 1){
		chomp($line);
		my @explode=split(/[$sep]/,$line);
		$string=$explode[$ensg]."_".$explode[$coord];
	
		$explode[$genotyper]=lc($explode[$genotyper]);
		$explode[$ensg]=lc($explode[$ensg]);
	
		if(!$identify{$string}){
			if($explode[$ensg]){
				if(!$hash{$explode[$ensg]}){
					$hash{$explode[$ensg]}=1;
					$hashDescription{$explode[$ensg]}=$explode[$description];#DESCRIPTION OF GENE
					$hashName{$explode[$ensg]}=$explode[$geneName];
					if($explode[$genotyper]=~/^P_Homo_var/){
						$hashGenotyper{$explode[$ensg]}=1;
					}
				} else {
					$hash{$explode[$ensg]}=$hash{$explode[$ensg]}+1;
					if($explode[$genotyper]=~/^P_Homo_var/) {
						$hashGenotyper{$explode[$ensg]}=1;
					}
				}
			}
			$identify{$string}=1;
		}
		if($explode[$ensg]){
			if(!$hash2{$explode[$ensg]}){
					$hash2{$explode[$ensg]}=1;
					$hashLines{$explode[$ensg]}[$hash2{$explode[$ensg]}]=$line;
					$hashScalar{$explode[$ensg]}=$hash2{$explode[$ensg]};
			}
		}
	}
}

foreach my $keys (keys %hash){
	my $keys2=uc($keys);
	print $DOM_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";

	for(my $i=0;$i<($hashScalar{$keys}+1) ;$i++){
		if($hashLines{$keys}[$i]){
			if($hashLines{$keys}[$i]=~/\W/){
				print $DOM_LINES "$hashLines{$keys}[$i]\n";
			}
		}
	}
	
	if($hash{$keys}>1){
		for(my $i=0;$i<($hashScalar{$keys}+1) ;$i++){
			if($hashLines{$keys}[$i]){
				if($hashLines{$keys}[$i]=~/\W/){
					print $REC_LINES "$hashLines{$keys}[$i]\n";
				}
			}
		}
		print $REC_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";
	} else {
		if($hashGenotyper{$keys}){
			if($hashGenotyper{$keys}==1){
				for(my $i=0;$i<($hashScalar{$keys}+1);$i++){
					if($hashLines{$keys}[$i]){
						if($hashLines{$keys}[$i]=~/\W/){
							print $REC_LINES "$hashLines{$keys}[$i]\n";
						}
					}
				}
				print $REC_HEADER "$keys2\t$hashName{$keys}\t$hash{$keys}\t$hashDescription{$keys}\n";
			}
		}
	}
}