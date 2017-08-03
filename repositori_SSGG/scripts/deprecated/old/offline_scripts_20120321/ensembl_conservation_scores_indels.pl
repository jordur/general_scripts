#!/usr/bin/perl -w

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

#
# Simple example to show how to get conservation scores for a slice.
# Works for ensembl release 51
#


sub usage
{
        print "\nEste script combina la información de las anotaciones del script SNP_effect.pl y la información original de las SNVs detectadas (output de samtools). Este script sólo añade las columnas del coverage y la calidad del SNP.\n";
        print "\nCOMO SE USA: sg_combinar_output_samtools_anotacion_SNP_effect.pl <fichero_samtools> <fichero_SNP_effect>\n";
        print "ejemplo: sg_combinar_output_samtools_anotacion_SNP_effect.pl chr10SNVs_hg19 SNVs_out_SNP_effect.txt\n\n";
        print "INPUT: Dos ficheros: el primero el fichero del samtools. El segundo fichero es el resultante de la anotación con SNP_effect.pl.\n\n";
        print "OUTPUT: fichero tabulado SNVs_descritas_ensembl_combinando_samtools_y_script_SNP_effect.txt con las siguientes columnas:Chr\tCoordinate\tReference\tGenotype\tGene_ID\tTranscript_ID\tConsequence\tPosition in cDNA\tSubstitution\tSNP_ID\tQ_consensus\tQ_SNV\tQ_max_map\tCoverage\n\n";
        exit(1);
}

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
        usage();
}

my @split_lineas;
my $contador=0;
open(INDELS,"<",$ARGV[0]);
while (my $lineas=<INDELS>)
{
        chomp($lineas);

	@split_lineas = split ("\t", $lineas);
#print $split_lineas[0],"\t",$split_lineas[1],"\t",$split_lineas[2],"\n";
my $reg = "Bio::EnsEMBL::Registry";
my $species = "Homo sapiens";
my $seq_region = $split_lineas[0];
my $dep;
if($split_lineas[1] - $split_lineas[2] == 0)
{
	$dep = 1;
}
else
{
	$dep = 0;
}
my $seq_region_start = $split_lineas[1];
my $seq_region_end =   $split_lineas[2] + $dep;
my $version = 59;
print  $split_lineas[0],"\t",$split_lineas[1],"\t", $split_lineas[2],"\t";
$reg->load_registry_from_db(
      -host => "ensembldb.ensembl.org",
      -user => "anonymous",
      -db_version => $version);

#get method_link_species_set adaptor
my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");

#get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");

throw("Unable to find method_link_species_set") if (!defined($mlss));

#get slice adaptor for $species
my $slice_adaptor = $reg->get_adaptor($species, 'core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);

#create slice 
my $slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region, $seq_region_start, $seq_region_end);
throw("No Slice can be created with coordinates $seq_region:$seq_region_start-$seq_region_end") if (!$slice);

#get conservation score adaptor
my $cs_adaptor = $reg->get_adaptor("Multi", 'compara', 'ConservationScore');			
#To get one score per base in the slice, must set display_size to the size of
#the slice.
my $display_size = $slice->end - $slice->start + 1; 
my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $display_size);

#print "number of scores " . @$scores . "\n";

#print out the position, observed, expected and difference scores.
my $resta=$split_lineas[2]-$split_lineas[1]+1;
#print "resto---- ",$resta,"\n";
foreach my $score (@$scores) {
    #$contador=$contador+1;
	if (defined $score->diff_score) 
	{
		if($dep == 1)
		{
#			print "entra en el primer if...score_position...",$score->position,"split_linea...",$split_lineas[1],"\n";
			if($score->position eq 1)
                       	{
				print $score->diff_score,"\n";
                        }
#			else
#			{	
#				print ",";
#			}
		}
		elsif($dep == 0)
		{
		print printf("%.3f", $score->diff_score),";";
		$contador=$contador+1;
		}
	}	
	#print "resto---- ",$resta,"\n";	
	elsif (defined $contador < defined $resta)
	{
		print "-------------------------------------------------";
	}
	
#print "--contador$contador--resta$resta\n"
}
print "\n";
$contador=0;
}
