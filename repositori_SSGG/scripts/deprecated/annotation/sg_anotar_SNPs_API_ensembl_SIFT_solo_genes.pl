#! /usr/bin/perl

#Descripción
#Script escrito en perl para realizar la anotación de los SNPs detectados en resecuenciaciones y análisis del exoma humano. Utiliza tres secciones del API de Ensembl: core, variation y compara. 
#Ubicación actual (cluster off-line):

#Consideraciones generales importantes
#-	Hay que realizar el análisis cromosoma a cromosoma.
#-	Este script funciona con la versión 37 del genoma humano (GRCh37 ó hg19).
#-	Cada línea corresponde a cada transcrito donde se encuentra el SNP. Así, podremos tener varias líneas por coordenada.

#En todos los archivos, el contenido de la columna viene determinado por la cabecera del mismo. Significados:
#chr: nombre del cromosoma.
#position:  la posición en el cromosoma.
#reference: base en la referencia.
#genotype: base detectada.
#Q_consensus: calidad del consenso.
#Q_SNV: calidad del SNP.
#Q_max_map: calidad máxima de mapeo.
#coverage: coverage.
#SNP_id: ID del SNP en la base de datos dbSNP del NCBI. Sólo se imprime si está descrito el SNP.
#gene: nombre del gen en el que se ha detectado el SNP.
#gene_description: descripción del gen en el que se ha detectado el SNP. No se imprime en todos los casos.
#strand:cadena del gen
#conservation: Este valor lo extraemos del API compara de Ensembl. Se basa en el valor GERP (http://genome.cshlp.org/content/15/7/901 ) de cada posición del genoma. En este script se imprime la diferencia entre el valor del cambio esperado menos el cambio observado. Esto implica que cuanto más negativo sea el valor impreso, el SNP tendrá más influencia sobre la secuencia. 

use strict;
use warnings;

use lib "/share/apps/src/ensembl-variation/modules/";
use lib "/share/apps/src/ensembl/modules/";
use lib "/share/apps/src/ensembl-compara/modules/";

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub usage
{
	print "\nCOMO SE USA: sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl <input> <cromosoma>\n\n";
	print "EJEMPLO: sg_anotar_SNPs_API_ensembl_SIFT_solo_genes.pl SNPs.txt 1\n\n";
	print "El argumento <cromosoma> debe ser un valor numérico. El formato 'chr1' no será reconocido\n\n";
	print "INPUT: El fichero de input debe ser un fichero tabulado con las siguientes columnas:\ncolumna 1: nombre del cromosoma (Sólo se admite el formato 'chr1')\ncolumna 2:  la posición en el cromosoma.\ncolumna 3: base en la referencia.\ncolumna 4: base detectada.\ncolumna 5: calidad del consenso.\ncolumna 6: calidad del SNP.\ncolumna 7: calidad máxima de mapeo.\ncolumna 8: coverage.\n\n";
	print "OUTPUT: obtenemos cuatro archivos:\n \nsnp_descritos_chr1.txt: listado de los SNPs que ya están descritos en NCBI/dbSNP.\nSIFT_snp_descritos_chr1.txt: listado de los SNPs descritos en formato para analizar por el programa SIFT (http://sift.jcvi.org/www/SIFT_chr_coords_submit.html).\nsnp_descritos_chr1.txt: listado de los SNPs que no están descritos en NCBI/dbSNP.\nSIFT_snp_no_descritos_chr1.txt: listado de los SNPs no descritos en formato para analizar por el programa SIFT.\n\n";
	exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}

# Decidimos con qué cromosoma vamos a trabajar
my $chr = $ARGV[1];
chomp($chr);

print "Intentando conectar con el servidor...";

# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
my $vfa = $reg->get_adaptor('human', 'variation', 'variationfeature');
my $sa = $reg->get_adaptor('human', 'core', 'slice');
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor("human", "core", "gene");

# Código de compara
#get method_link_species_set adaptor
my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
#get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE","mammals");

my $slice_adaptor_compara = $reg->get_adaptor("Homo sapiens",'core', 'Slice');
#get conservation score adaptor
my $cs_adaptor = $reg->get_adaptor("Multi", 'compara', 'ConservationScore');

print "\n\nConectado!\n";

print "\nAnalizando los datos...\n";

# connect to Variation database
my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
  -host    => 'ensembldb.ensembl.org',
  -user    => 'anonymous',
  -species => 'Homo sapiens',
  -group   => 'variation',
  -dbname =>  $vfa->dbc->dbname,
  -port => '5306'
);

# connect to Core database

my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host    => 'ensembldb.ensembl.org',
  -user    => 'anonymous',
  -species => 'Homo sapiens',
  -group   => 'core',
  -dbname =>  $sa->dbc->dbname,
  -port => '5306'
);

print "dbname_variation:",$vfa->dbc->dbname,".......dbname_core.....",$sa->dbc->dbname,"\n";

# Abrimos los archivos de salida
#
# Abrimos el archivo de salida donde se recogen todos los SNPs ya descritos
open(DESCRITOS,'>',"snp_descritos_chr$chr.txt");
print DESCRITOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","gene_name","\t","gene_descrition","\t","strand","\t","conservation","\n";

# Abrimos el archivo de salida donde se recogen todos los SNPs ya descritos en el formato para SIFT
open(DESCRITOS_SIFT,'>',"SIFT_snp_descritos_chr$chr.txt");

# Abrimos el archivo de salida donde se recogen todos los SNPs no descritos
open(NO_DESCRITOS,'>',"snp_no_descritos_chr$chr.txt");
print NO_DESCRITOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","gene_name","\t","gene_descrition","\t","strand","\t","conservation","\n";

# Abrimos el archivo de salida donde se recogen todos los SNPs no descritos en el formato para SIFT
open(NO_DESCRITOS_SIFT,'>',"SIFT_snp_no_descritos_chr$chr.txt");

# Abrimos el archivo de los SNPs

open(SNP,'<',$ARGV[0]);

while (my $snp=<SNP>)
{
        chomp($snp);
        my @ex = split(/\t/,$snp);
	my $chr_input = "chr".$chr;
	
	if($chr_input eq $ex[0])
	{
		my $start_slice = $ex[1];
		my $end_slice = $ex[1];
		my $slice_adaptor = $dbCore->get_SliceAdaptor(); #get the database adaptor for Slice objects
		my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$start_slice,$end_slice);
#		my $slice_reverse = $slice_adaptor->fetch_by_region('chromosome',$chr,$start_slice,$end_slice,'-1');
		my $slice_chr = $slice_adaptor->fetch_by_region('chromosome',$chr);
		# create a new VariationFeature object
		my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor();
                my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
                -start => $ex[1],
                -end => $ex[1],
                -slice => $slice_chr,           # the variation must be attached to a slice
                -allele_string => $ex[2],"/",$ex[3],    # the first allele should be the reference allele
#                -strand => 1,
                -map_weight => 1,
                -source  => 'dbSNP',
                -adaptor => $vf_adaptor,           # we must attach a variation feature adaptor
                );
		
		# Comprobamos si el SNP está descrito
		my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice
                my $rsID = codigo_NCBI($vfs);
#		print "rsID____$rsID\n";
		
		# Si el SNP está descrito, lo imprimimos en el fichero de los descritos, sinó, lo imprimimos en el de los no descritos

		# Convertimos la base del genotipo a homocigota
		my $base_homo = convertir_a_homocigoto($ex[2],$ex[3]);

		if($rsID ne "N/A")
		{
			# Calculamos el score de GERP con la subrutina:

                        my $GERP = GERP_score($ex[1],$chr,$mlss);

			foreach my $gene(@{$gene_adaptor->fetch_all_by_Slice($slice)})
                        {
                                print DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\t",$gene->external_name,"\t",$gene->description,"\t",$gene->strand(),"\t",$GERP,"\n";
				print DESCRITOS_SIFT "$chr",",",$ex[1],",",$gene->strand(),",",$ex[2],"/",$base_homo,"\n";
                        }

#			print DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t","\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\n";
		}
		else
		{
			# Calculamos el score de GERP con la subrutina:
			
			my $GERP = GERP_score($ex[1],$chr,$mlss);

			foreach my $gene(@{$gene_adaptor->fetch_all_by_Slice($slice)})
	                {
				
        	                print NO_DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\t",$gene->external_name,"\t",$gene->description,"\t",$gene->strand(),"\t",$GERP,"\n";
				print NO_DESCRITOS_SIFT "$chr",",",$ex[1],",",$gene->strand(),",",$ex[2],"/",$base_homo,"\n";
                	}

#	                foreach my $gene_reverse(@{$gene_adaptor->fetch_all_by_Slice($slice_reverse)})
#        	        {
#                	         print NO_DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\t",$gene_reverse->external_name,"\t",$gene_reverse->description,"\t",$slice_reverse->strand(),"\t",$GERP,"\n";
#                	}
		}
	}
}


# Cerramos el archivo con todos los SNPs
close(DESCRITOS);
close(DESCRITOS_SIFT);
close(NO_DESCRITOS);
close(NO_DESCRITOS_SIFT);
close(SNP);


# Salimos del programa
exit;


#####################################################################################################################


#						SUBRUTINAS


####################################################
#
# Subrutina para extraer el número de cromosoma
#
#	NO LA UTILIZAMOS PORQUE ESTÁ MAL
#
####################################################

sub def_chr
{
        if(length(@_) > 3)
        {
		print "@....@_\n";
                return substr(@_,3,length(@_)-3);
        }
        else
        {
		print "@....@_\n";
                return @_;
        }
}

#####################################################


#####################################################
#
# Subrutina para extraer el código NCBI del SNP
#
#####################################################

sub codigo_NCBI
{
        my $no = 'N/A';
        foreach my $vf (@{$_[0]})
        {
                my $codigo = $vf->variation_name;
                if(length($codigo) > 4)
                {
                        return $codigo;
                }
        }
        return $no;
}


#####################################################

######################################################################################################
#
#               Subrutina de compara para obtener el valor de GERP
#
######################################################################################################

sub GERP_score
{
        my $no_GERP = 'unknown';
        # creando el slice para calcular el valor GERP
        my $intervalo = 1;
        my $start_compara = $_[0] - $intervalo;
        my $end_compara = $_[0] + $intervalo;
        my $slice_compara = $slice_adaptor_compara->fetch_by_region('toplevel', $_[1], $start_compara, $end_compara);
        #To get one score per base in the slice, must set display_size to the size of the slice.
        my $display_size = $slice_compara->end - $slice_compara->start + 1;
        my $medio = $intervalo + 1;
        my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($_[2], $slice_compara, $display_size);

        #print out the position, observed, expected and difference scores.
        foreach my $score (@$scores)
        {
                if (defined $score->diff_score)
                {
                        if($score->position eq $medio)
                        {
                                return my $return = $score->diff_score;
                        }
                }
        }
        return $no_GERP;
}

################################################################################################


##########################################################################################
#
# Subrutina convertir a homocigoto
#
##########################################################################################

# Necesita dos variables. La primera debe ser la referencia y la segunda el genotipo. No acepta IUPAC V, H, D, B, X o N

sub convertir_a_homocigoto
{
        my $out;
        if(($_[1] eq "A") || ($_[1] eq "C") || ($_[1] eq "G") || ($_[1] eq "T"))
        {
                return $_[1];
        }
        else
        {
                if($_[1] eq "R")
                {
                        if($_[0] eq "A")
                        {
                                $out = "G";
			}
                        else
                        {
                                $out = "A";
                        }
                        return $out;
                }
                elsif($_[1] eq "M")
                {
                        if($_[0] eq "A")
                        {
                                $out = "C";
                        }
                        else
                        {
                                $out = "A";
                        }
                        return $out;
                }
                elsif($_[1] eq "W")
                {
                        if($_[0] eq "A")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "A";
                        }
                        return $out;
                }
                elsif($_[1] eq "S")
                {
                        if($_[0] eq "C")
                        {
                                $out = "G";
                        }
                        else
                        {
                                $out = "C";
                        }
                        return $out;
                }
                elsif($_[1] eq "Y")
		{
                        if($_[0] eq "C")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "C";
                        }
                        return $out;
                }
                elsif($_[1] eq "K")
                {
                        if($_[0] eq "G")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "G";
                        }
                        return $out;
                }
                else
                {
                die("este script no acepta IUPAC V, H, D, B, X o N");
                }
        }
}

