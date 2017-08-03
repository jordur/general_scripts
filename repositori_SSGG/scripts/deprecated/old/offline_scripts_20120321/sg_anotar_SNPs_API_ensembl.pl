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
#consequence: consecuencia del SNP.
#transcript_ID: ID de Ensembl del transcrito.
#SNP_id: ID del SNP en la base de datos dbSNP del NCBI. Sólo se imprime si está descrito el SNP.
#gene: nombre del gen en el que se ha detectado el SNP.
#gene_description: descripción del gen en el que se ha detectado el SNP. No se imprime en todos los casos.
#HGVS: nomenclatura del SNP en formato HGVS (http://www.hgvs.org/mutnomen/recs-DNA.html ). La estructura es la siguiente:
#ID del trancritog.positiongenotype>reference 
#conservation: Este valor lo extraemos del API compara de Ensembl. Se basa en el valor GERP (http://genome.cshlp.org/content/15/7/901 ) de cada posición del genoma. En este script se imprime la diferencia entre el valor del cambio esperado menos el cambio observado. Esto implica que cuanto más negativo sea el valor impreso, el SNP tendrá más influencia sobre la secuencia. 
#Q_consensus: calidad del consenso.
#Q_SNV: calidad del SNP.
#Q_max_map: calidad máxima de mapeo.
#coverage: coverage.

use strict;
use warnings;

use lib "/share/apps/src/ensembl-variation/modules/";
use lib "/share/apps/src/ensembl/modules/";
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub usage
{
	print "\nCOMO SE USA: sg_anotar_SNPs_API_ensembl.pl <input> <cromosoma>\n\n";
	print "EJEMPLO: sg_anotar_SNPs_API_ensembl.pl SNPs.txt 1\n\n";
	print "El argumento <cromosoma> debe ser un valor numérico. El formato 'chr1' no será reconocido\n\n";
	print "INPUT: El fichero de input debe ser un fichero tabulado con las siguientes columnas:\ncolumna 1: nombre del cromosoma (Sólo se admite el formato 'chr1')\ncolumna 2:  la posición en el cromosoma.\ncolumna 3: base en la referencia.\ncolumna 4: base detectada.\ncolumna 5: calidad del consenso.\ncolumna 6: calidad del SNP.\ncolumna 7: calidad máxima de mapeo.\ncolumna 8: coverage.\n\n";
	print "OUTPUT: obtenemos siete archivos:\n snp_ensembl_chr1.txt: listado de la anotación de todos los SNPs.\nsnp_descritos_chr1.txt: listado de los SNPs que ya están descritos en NCBI/dbSNP.\nsnp_intergenic_chr1.txt: listado de los SNPs en regiones intergénicas.\nsnp_no_sinonimos_chr1.txt: listado de los SNPs cuya consecuencia es un cambio no sinónimo.\nsnp_splicing_chr1.txt: listado de los SNPs que pueden afectar al splicing.\nsnp_no_codificantes_chr1.txt: listado de los SNPs en regiones no codificantes.\nsnp_sinonimos_chr1.txt: listado de los SNPs cuya consecuencia es un cambio sinónimo.\n\n";
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
my $exon_adaptor = $reg->get_adaptor( 'Human', 'Core', 'Exon' );

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
#  -dbname  => 'homo_sapiens_variation_58_37c'
  -port => '5306'
);

# connect to Core database

my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host    => 'ensembldb.ensembl.org',
  -user    => 'anonymous',
  -species => 'Homo sapiens',
  -group   => 'core',
  -dbname =>  $sa->dbc->dbname,
#  -dbname  => 'homo_sapiens_core_58_37c'
  -port => '5306'
);

# Abrimos los archivos de salida
#
# Abrimos el archivo de salida donde se recogen todos los SNPs ya descritos
open(DESCRITOS,'>',"snp_descritos_chr$chr.txt");
print DESCRITOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\n";

# Abrimos el archivo de salida donde se recogen todos los SNPs no descritos
open(NO_DESCRITOS,'>',"snp_no_descritos_chr$chr.txt");
print NO_DESCRITOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";


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
		my $slice_chr = $slice_adaptor->fetch_by_region('chromosome',$chr);
		# create a new VariationFeature object
		my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor();
                my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
                -start => $ex[1],
                -end => $ex[1],
                -slice => $slice_chr,           # the variation must be attached to a slice
                -allele_string => $ex[2]."/".$ex[3],    # the first allele should be the reference allele
                -strand => $ex[4],
                -map_weight => 1,
                -source  => 'dbSNP',
                -adaptor => $vfa,           # we must attach a variation feature adaptor
                );
		
		# Comprobamos si el SNP está descrito
		my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice
                my $rsID = codigo_NCBI($vfs);
#		print "rsID____$rsID\n";
		
		# Si el SNP está descrito, lo imprimimos en el fichero de los descritos, sinó, lo imprimimos en el de los no descritos
	
		if($rsID ne "N/A")
		{
			print DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t","N/A","\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\n";
		}
		else
		{
			foreach my $con(@{$new_vf->get_all_TranscriptVariations})
                	{
				# Calculamos el score de GERP con la subrutina:
			
				my $GERP = GERP_score($ex[1],$chr,$mlss);

				foreach my $string(@{$con->consequence_type})
                        	{
					print "coordenada...$ex[1]......$string\n";
                                	# Imprimimos en un fichero los SNPs que están en zonas intergénicas
                                	if ($string eq "INTERGENIC")
                                	{
						print NO_DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t",$string,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\t---\t---\t---\tg.",$ex[1],$ex[3],">",$ex[2],"\t",$GERP,"\n";
                                	}
                                	else
                                	{
					my $transcript = $con->transcript->stable_id;
#        		                print "stable id..... $transcript\n";

                        		# Obtenemos el nombre del gen
                                	my $gene = $gene_adaptor->fetch_by_transcript_stable_id($transcript);
#                                	print "gene..... $gene\n";
                                	my $external_name = $gene->external_name;
                                	my $description = $gene->description();

                                	# Coordenadas del cDNA y de la zona codificante-CDS
                                	my $cdna_start = $con->transcript->cdna_coding_start;
                                	my $cdna_end = $con->transcript->coding_region_end;
                                	# Coordenadas de la variación
                                	my $coord_cDNA = $con->_calc_transcript_coords;
                                        print NO_DESCRITOS "chr$chr\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$rsID,"\t",$string,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\t",$transcript,"\t",$external_name,"\t",$description,"\t",$transcript,"g.",$ex[1],$ex[3],">",$ex[2],"\t",$GERP,"\t","\n";
					}
				}
			}
		}
	}
}


# Cerramos el archivo con todos los SNPs
close(DESCRITOS);
close(NO_DESCRITOS);
close(SNP);


# INTERGÉNICOS: Abrimos el archivo de salida donde se recogen los SNPs que están zonas intergénicas
open(INTERGENIC,">","snp_intergenic_chr$chr.txt");
print INTERGENIC "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";

# SINÓNIMOS: Abrimos el archivo de salida donde se recogen los SNPs no descritos en el NCBI y cuyos cambios son sinónimos
open(SINONIMOS,">","snp_sinonimos_chr$chr.txt");
print SINONIMOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";

# NO SINÓNIMOS: Abrimos el archivo de salida donde se recogen los SNPs no descritos en el NCBI y cuyos cambios no son sinónimos
open(NOSINONIMOS,">","snp_no_sinonimos_chr$chr.txt");
print NOSINONIMOS "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";

# NO CODIFICANTES: Abrimos el archivo de salida donde se recogen los SNPs no descritos en el NCBI y que se encuentran en zonas no codificantes
open(NOCODIF,">","snp_no_codificantes_chr$chr.txt");

# SPLICING: Abrimos el archivo de salida donde se recogen los SNPs no descritos en el NCBI y que se encuentran en zonas de splicing
open(SPLICING,">","snp_splicing_chr$chr.txt");
print SPLICING "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";

# UTR: Abrimos el archivo de salida donde se recogen los SNPs no descritos en el NCBI y que se encuentran en zonas 5' y 3' UTR
open(UTR,">","snp_UTR_chr$chr.txt");
print UTR "chr","\t","coordinate","\t","reference","\t","genotype","\t","SNP_ID","\t","consequence","\t","Q_consensus","\t","Q_SNV","\t","Q_max_map","\t","coverage","\t","transcript_ID","\t","gene_name","\t","gene_descrition","\t","SNV_annotation(HGVS)","\t","conservation","\n";

# Imprimimos los ficheros organizando los SNPs por consecuencia

open(TODO,'<',"snp_no_descritos_chr$chr.txt");

while (my $anotados=<TODO>)
{
        chomp($anotados);
        my @snps_anotados = split(/\t/,$anotados);
        if($snps_anotados[4] eq "INTERGENIC")
        {
                print INTERGENIC $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
        }
        else
        {
        	if($snps_anotados[4] eq "SYNONYMOUS_CODING")
        	{
			print SINONIMOS $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
		}
		elsif($snps_anotados[4] eq "NON_SYNONYMOUS_CODING")
                {
			print NOSINONIMOS $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
                }
                elsif(($snps_anotados[4] eq "ESSENTIAL_SPLICE_SITE") || ($snps_anotados[4] eq "SPLICE_SITE"))
                {
			print SPLICING $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
                }
		elsif(($snps_anotados[4] eq "5PRIME_UTR") || ($snps_anotados[4] eq "3PRIME_UTR"))
		{
			print UTR $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
		}
                else
                {
			print NOCODIF $snps_anotados[0],"\t",$snps_anotados[1],"\t",$snps_anotados[2],"\t",$snps_anotados[3],"\t",$snps_anotados[4],"\t",$snps_anotados[5],"\t",$snps_anotados[6],"\t",$snps_anotados[7],"\t",$snps_anotados[8],"\t",$snps_anotados[9],"\t",$snps_anotados[10],"\t",$snps_anotados[11],"\t",$snps_anotados[12],"\t",$snps_anotados[13],"\t",$snps_anotados[14],"\n";
                }
        }
}

close(INTERGENIC);
close(SINONIMOS);
close(NOSINONIMOS);
close(NOCODIF);
close(SPLICING);
close(UTR);
close(TODO);



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




######################################################################################################
#
#               Subrutina para obtener el código del transcrito
#
#               NO LA GASTAMOS!!
#
######################################################################################################

sub transcript_con
{
        my $no_transcript = 'no_gene';
        #print out the position, observed, expected and difference scores.
        my $trans = $_[0]->transcript->stable_id;
                if(length($trans) > 7)
                {
                        return $trans;
                }
        return $no_transcript;
}

##########################################################################################################
