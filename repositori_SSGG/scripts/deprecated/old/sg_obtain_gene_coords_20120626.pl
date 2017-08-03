#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for obtaining gene coordinates from th HGNC names
# @Author: Arbol
# @Contributors: biomart & ensembl scripts
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use File::Basename;
use Cwd;
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;
use warnings;
use Text::CSV;

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for obtaining gene coordinates from the Ensembl-biomart databases.

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 input		File containing the gene names in HGNC notation and the desired transcript names for the given gene (in ENSEMbl notation), in CSV format (comma separated values, withou spaces).
 			If no transcripts are included for a given gene, then all the transcripts exons will be output.
 			Example of contents of input file:
 				ABCC9,ENST00000544039,ENST00000261201
				ACTC1,ENST00000290378,ENST00000557860,ENST00000544062
				ACTN2

* Options: *
 -prefix sampleprefix_ 	Adds the given string as prefix for the output results
 -outdir output_path	Path where the results will be output. Default path is <cwd>/panel_info, where <cwd> is the current working directory

* Examples of running: *
 \$$scriptname -outdir arbol_mola -prefix BM42_ panel_CSV.txt
 \n";

	exit(1);
}

sub add_and_short_intervals {
	# This subroutine shorts and merges new intervals into an existent intervals array.
	
	# Arguments:
	# Please note that $new_interval, $intervals_array and $added_pos are all references!!
	my ( $intervals_array, $new_interval )=@_;
	
	# Normalize current interval (low,high):
	my @current_interval;
	if ( ${$new_interval}[1] > ${$new_interval}[0] ){ @current_interval = (${$new_interval}[0], ${$new_interval}[1]); }
	else { @current_interval = (${$new_interval}[1], ${$new_interval}[0]); }
							
	# The exon coordinates intervals get ordered in the @exon_intervals array.
	# These intervals get also merged with adjacent intervals
						
	# The first interval is directly added to the array
	if ( $#{$intervals_array} == -1 ) {
		push (@{$intervals_array}, \@current_interval);
	}
	else {
		# Check if the exon interval intersects with existent intervals,
		# and in which position it should be added
		my $low_pos = 0;
		while ( $low_pos <= $#{$intervals_array} && $current_interval[0] >= ${$intervals_array}[$low_pos][0] ) { $low_pos ++; }
		if ( $low_pos == 0 ) {
			# Add new interval at the beginning of the intervals_array
			if ( $current_interval[1] < ${$intervals_array}[0][0] ) {
				unshift(@{$intervals_array}, \@current_interval);
			}
			# Modify the low_pos of intevals_array
			else {
				${$intervals_array}[0][0]=$current_interval[0];
			}
		} else {
			# Modify the high_pos of intervals_array
			if ( $current_interval[0] <= ${$intervals_array}[$low_pos-1][1] ) {
				${$intervals_array}[$low_pos-1][1]=$current_interval[1]; 
			}
			# Add new interval at $low_pos position of the intervals_array
			else { 
				splice(@{$intervals_array}, $low_pos, 0, \@current_interval); 
			}
		}
	}
}

sub gene_coords_parser {
	# This subroutine parses the contents of the output file created by the Biomart query, 
	# checks genes codes & names, duplicates, etc.
	# Finally, it writes the results to a .txt file (tab separated values) that can be uploaded
	# to Agilent's eArray
	# Another file containing also the splicing zones (30nt at each side of the interval) is created
	
	# Arguments:
	# $genes_file: csv file containing all the exon coordinates (for all transcripts of the gene)
	# $gene_name: name of the gene to consider
	# $output_file: name of the file to output the normalized exon coordinates for desired transcripts
	# $transcript: list of the transcripts to consider from gene. All transcripts are considered if empty. This parameter is passed as a reference
	# $splice_file: name of the file to output the normalized exon coordinates including splicing sites (30 nt at each side of the exon)
	# $sum_exons: variable to return the amount of bases covered by the normalized exon intervals
	# $sum_splice_exons: variable to return the amount of fases covered by the normalized exon intervals including splicing sites
	my ($genes_file,$gene_name,$output_file,$transcript,$splice_file,$sum_exons,$sum_splice_exons) = @_;
		
	# variables:
	my $csv = Text::CSV->new();
	my @biomart_row;
	my $row = 0;
	my @exon_intervals; # Array containing shorted & merged exon coordinates intervals
	my @splice_exon_intervals; # Array containing shorted & merged exon coordinates intervals
	my $chromosome = "";
	
	open (CSV, "<", $genes_file) or die("Couldn't open file $genes_file!! $!\n");
	
	while (<CSV>) {
		$row ++;
		if ( $row > 1 ) {
			if ($csv->parse($_)) {
				@biomart_row = $csv->fields();
				# Check that the gene name and chromosome is the same for all exons
				if ($chromosome eq "") { $chromosome = $biomart_row[6]; }
				if ($biomart_row[7] ne $gene_name || $biomart_row[6] ne $chromosome ) {
					print "Problem in row $row del fichero $genes_file. Gene name ($biomart_row[8]) not expected (different to $gene_name).\n";			
				} else {
					# Check all the transcripts and include their exon coordinates in the output file
					my $trans_amount=scalar(@{$transcript});
					
					# Variable to store the added positions from the new interval added to the exon_intervals array
					my $added_coords=0;
					
					# If NO transcript has been defined, then all exon coordinates will be considered.
					# Otherwise, only the exon coordinates from the desired transcripts
					if ($trans_amount == 0) {
						my @new_interval=($biomart_row[8], $biomart_row[9]);
						my @new_splice_interval;
						if ( $new_interval[0] > 30 ) { $new_splice_interval[0] = $new_interval[0] - 30; }
						else { $new_splice_interval[0] = 0; }
						$new_splice_interval[1] = $new_interval[1] + 30;
						add_and_short_intervals( \@exon_intervals, \@new_interval );
						add_and_short_intervals( \@splice_exon_intervals, \@new_splice_interval );
					} else {
						for (my $i=0; $i < $trans_amount; $i++) {
							if ( $biomart_row[1] eq ${$transcript}[$i] ) {
								my @new_interval=($biomart_row[8], $biomart_row[9]);
								my @new_splice_interval;
								if ( $new_interval[0] > 30 ) { $new_splice_interval[0] = $new_interval[0] - 30; }
								else { $new_splice_interval[0] = 0; }
								$new_splice_interval[1] = $new_interval[1] + 30;								
								add_and_short_intervals( \@exon_intervals, \@new_interval );
								add_and_short_intervals( \@splice_exon_intervals, \@new_splice_interval );
							}
						}
					}		
				}
			} else {
				my $err = $csv->error_input;
				print "Failed to parse line $row: $err\n";
			}
		}
	}
	
	close CSV;
	open (OUTPUT, ">", $output_file) or die("Couldn't open output file $output_file!! $!\n");
	open (SPLICEFILE, ">", $splice_file) or die("Couldn't open output file $splice_file!! $!\n");
	
	# Normalization of all intervals for eArray (Agilent's utility for designing panels)
	# Exon intervals under 60 nt will be "enlarged" until reaching this size
	# The splicing zones are calculated adding 30nt at each side of the exon intervals
	
	$$sum_exons = 0;
	$$sum_splice_exons = 0;
		
	my $sum_desired_exons = 0;
	# Chromosome name "MT" is renamed to "M"
	if ( $chromosome eq "MT" ) { $chromosome = "M"; }
	for (my $i=0; $i < $#exon_intervals; $i++) {
		# Splicing zones added to intervals:
		my $out_row = "chr" . $chromosome . ":" . $splice_exon_intervals[$i][0] . "-" . $splice_exon_intervals[$i][1] . "\n";
		print SPLICEFILE $out_row;
		$$sum_splice_exons = $$sum_splice_exons + $splice_exon_intervals[$i][1] - $splice_exon_intervals[$i][0];

		# Exon intervals:
		$out_row = "chr" . $chromosome . ":" . $exon_intervals[$i][0] . "-" . $exon_intervals[$i][1] . "\n";
		print OUTPUT $out_row;
		$$sum_exons = $$sum_exons + $exon_intervals[$i][1] - $exon_intervals[$i][0];
	}
	
	close SPLICEFILE;
	close OUTPUT;
	
	# Display error if no gene with the given name was found!
	if ( $row < 2 ) {
		print "Problem with gene $gene_name!! No information found!!\n";
	}
	
	return $sum_desired_exons;
}

# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
my $fooGlobalVar = $ENV{"PATH"};

# ------------------------
# -- default parameters --
# ------------------------
my $prefix = "";
my $outdir = cwd() . "/panel_info";
my $genes_file = "";
my $transcript = "";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.1
Script for obtaining gene coordinates from the Ensembl-biomart databases · Copyright 2012 by Arbol 
**************************************************************************************************\n";

# If no arguments are passed to the script, the script usage gets displayed
if ($#ARGV+1 < 1) {
	usage( $basename.$ext );
}

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------
while ($#ARGV >= 0) {
	my $opt   = shift @ARGV;
	if ($opt =~ /^\-/) {
		if ($opt eq "-outdir") {
			$outdir = shift @ARGV;
        }
		elsif ($opt eq "-prefix") {
			$prefix = shift @ARGV;
        }
		else {
			print "Bad option: $opt\n";
			print "Call $basename$ext with no parameters for help.\n";
			exit(1);							
		}
	}
	else {
		# Mandatory parameters definitions
		$genes_file = $opt;
	}
}
if ($genes_file eq "") {
	print "Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

my $now = localtime;
print "Script started at $now\n";

# Parse the input CSV files with genes names and transcripts
my $csv = Text::CSV->new();
my @gene_transcripts;
my @genes_array;

open (CSV, "<", $genes_file) or die $!;
my $row = 0;
while (<CSV>) {
	if ($csv->parse($_)) {
		@gene_transcripts = $csv->fields();
		print "Gene: $gene_transcripts[0] Transcripts: @gene_transcripts[1..$#gene_transcripts]\n";
		push (@{$genes_array[$row]}, @gene_transcripts);
		$row ++;
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}
}
close CSV;

# Create output folder:
unless(-d $outdir){
	mkdir $outdir or die ("Impossible to create folder $outdir!!");
}

# Create and open stats file:
my $stats_file=$outdir . "/". $prefix . "Stats_exon.txt";
open STATS, '>', $stats_file or die ("Impossible to create $stats_file!! $!");

# Variables for the total number of bases of the design (without and including splice sites):
my $all_bases = 0;
my $all_bases_including_splice = 0;

# Auxiliar variables for storing all the files names (needed for concatenation of these files contents):
my @all_splice_files;
my @all_exons_files;

# Obtain coordinates and parameters for each gene
for ( my $i=0; $i<=$#genes_array; $i++ ) {
	
	# The whitespaces from the beginning of the gene names should be removed:
	my $gene_name = $genes_array[$i][0];
	$gene_name =~ s/^\s+//;
	my @gene = [$gene_name];
	
	my $output_file = $outdir . "/". $prefix . $gene_name . "_all_isoforms.csv";
	print "\nRetrieving results for gene $gene_name to file $output_file\n";
	
	# This perl API representation is only available for configuration versions >=  0.5 
	
	#"PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/.
	my $confFile = "/share/apps/src/biomart-perl/conf/ensembl.xml";  
	# For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry";
	#
	# NB: change action to 'clean' if you wish to start a fresh configuration  
	# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
	#
	
	my $action='cached';
	my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
	my $registry = $initializer->getRegistry;
	
	my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
		
		$query->setDataset("hsapiens_gene_ensembl");
		$query->addFilter("hgnc_symbol", @gene);
		$query->addAttribute("ensembl_gene_id");
		$query->addAttribute("ensembl_transcript_id");
		$query->addAttribute("5_utr_start");
		$query->addAttribute("5_utr_end");
		$query->addAttribute("3_utr_start");
		$query->addAttribute("3_utr_end");
		$query->addAttribute("chromosome_name");
		$query->addAttribute("external_gene_id");
		$query->addAttribute("exon_chrom_start");
		$query->addAttribute("exon_chrom_end");
	
	$query->formatter("CSV");
	
	
	my $query_runner = BioMart::QueryRunner->new();
	############################## GET COUNT ############################
	# $query->count(1);
	# $query_runner->execute($query);
	# print $query_runner->getCount();
	#####################################################################
	
	# Redirect STDOUT to output_file:
	{
		# Open the output file for the results (it will capture STDOUT)
		local *STDOUT;
		open(STDOUT, '>', $output_file) || die("Can't redirect STDOUT: $!");
		
		############################## GET RESULTS ##########################
		# to obtain unique rows only
		# $query_runner->uniqueRowsOnly(1);
		
		$query_runner->execute($query);
		$query_runner->printHeader();
		$query_runner->printResults();
		$query_runner->printFooter();
		#####################################################################

		# close the file for outputting the results:
		close STDOUT;
	}
	my $cols = @{ $genes_array[$i] };
	my @transcripts = @{ $genes_array[$i] } [ 1..($cols-1) ];
	my $exons_file = $outdir . "/". $prefix . $gene_name . "_exon_coords.txt";
	my $splice_file = $outdir . "/". $prefix . $gene_name . "_splice+exon_coords.txt";
	my $sum_desired_exons;
	my $sum_desired_splice_exons;
	# Please note that in the call to gene_coords_parser, some parameters are passed as references!!
	gene_coords_parser( $output_file, $gene_name, $exons_file, \@transcripts, $splice_file, \$sum_desired_exons, \$sum_desired_splice_exons );
	my $line = "Gene " . $gene_name .":\n";
	print STATS $line;
	$line = "Quantitty of bases from desired transcripts (" . @transcripts . "): " . $sum_desired_exons . " bp\n";
	print STATS $line;
	$line = "Quantitty of bases from desired transcripts including splice sites (" . @transcripts . "): " . $sum_desired_splice_exons . " bp\n";
	print STATS $line;
	
	$all_bases = $all_bases + $sum_desired_exons;
	$all_bases_including_splice = $all_bases_including_splice + $sum_desired_splice_exons; 
	
	# The exons and splice+exons files will be concatenated in a single file. For this reason they are added to an array, so that they can later
	# be processed after getting them in the ARGV array
	push (@all_splice_files, $splice_file);
	push (@all_exons_files, $exons_file);
}

my $line = "------------------------------------------------------------------- \n";
print STATS $line;
$line = "Total number of bases of the design: $all_bases bp \n";
print STATS $line;
$line = "Total number of bases of the design including splice sites: $all_bases_including_splice bp \n";
print STATS $line;
close STATS;

print "Total number of bases of the design: $all_bases bp \n";
print "Total number of bases of the design including splice sites: $all_bases_including_splice bp \n";

# Create a single file with all the exons files:
my $design_exons_file=$outdir . "/". $prefix . "All_exon_coords.txt";
@ARGV = @all_exons_files;
open SEL, '>', $design_exons_file or die ("Impossible to create $design_exons_file!! $!");

# Print all exons coords in the file
while (<>) {
	print SEL;
}
close SEL;

# Create a single file with all the splice + exons files:
my $design_splice_file=$outdir . "/". $prefix . "All_splice+exon_coords.txt";
@ARGV = @all_splice_files;
open SEL, '>', $design_splice_file or die ("Impossible to create $design_splice_file!! $!");

# Print all exons coords in the file
while (<>) {
	print SEL;
}
close SEL;

$now = localtime;
print "Script finished at $now\n";

#******************************************************************* 
