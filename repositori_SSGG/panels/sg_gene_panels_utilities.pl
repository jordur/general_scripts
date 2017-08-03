#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script for obtaining gene coordinates from HGNC names or Ensembl IDs. It also checks and modifys existent design panels
# @Author: Arbol
# @Contributors: biomart & ensembl scripts
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use File::Basename;
use Cwd;
use strict;
use warnings;
use Data::Dumper;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;
use Text::CSV;
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];
use POSIX qw( strftime );
use File::Copy;
use intervals qw( add_and_short_intervals subtract_intervals print_and_count_intervals add_intervals_to_hash add_intervals_hashes subtract_intervals_hashes normalize_intervals_hash obtain_splicing_from_exons );
use transcripts qw( new_transcript show_transcript );
use bioinfo_files_utilities qw( import_intervals_file create_intervals_file );

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for obtaining gene coordinates from the Ensembl-biomart databases, and to check its results against designs (Agilent, Illumina...), allowing to change the panel coordinates for submitting new designs

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 input		File containing the gene names in HGNC notation and the desired transcript names for the given gene (in ENSEMbl notation), in CSV format (comma separated values, without spaces).
 			If no transcripts are included for a given gene, then all the transcripts exons will be output.
 			Example of contents of input file:
 				ABCC9,ENST00000544039,ENST00000261201
				ACTC1,ENST00000290378,ENST00000557860,ENST00000544062
				ACTN2

* Options: *
** General options **
 -prefix sampleprefix_     Adds the given string as prefix for the output results
 -biomart conf_file        Path to the biomart config file (default /share/apps/src/biomart-perl/conf/ensembl.xml)

** For obtaining designs: **
 -outdir output_path	   Path where the results will be output (default path is <cwd>/panel_info, where <cwd> is the current working directory)
 -splicing #bases          Number of bases before and after the exon to consider as splicing sites (default 20)
 -minsize #bases           Minimum number of bases of the panels capture regions (default 0)
 -no_ensembl_UTR5          Do not include known Ensembl 5' UTRs (by default included)
 -no_ensembl_UTR3          Do not include known Ensembl 3' UTRs (by default included)
 -only_biobase_UTR5        Only include 5' UTRs from transcripts described in BioBase
 -only_biobase_UTR3        Only include 3' UTRs from transcripts described in BioBase
 -UTR5 #bases              Number of bases to consider for the 5' UTR in case that are not found in ensembl (default 158, 75% of known 5' UTRs threshold)
 -UTR3 #bases              Number of bases to consider for the 3' UTR in case that are not found in ensembl (default 502, 75% of known 3' UTRs threshold)
 -intervals tsv_file       Tab separated values (tsv, chr \t pos1 \t pos2) file containing other intervals to be included in the design (for example, gene expansion zones in introns, ...)
 -attempts #attempts       Maximum number of attempts when querying Ensembl-biomart (default 3) 

** For checking capture regions or panels (SureSelect, HaloPlex, TruSeq, ...): **
 -target design_folder	   Folder containing the target results of previous designs, which will be checked against the desired panel (SureSelect, HaloPlex, TruSeq, etc...)
 -check capture_file	   File containing the capture region (panel) coordinates (from Agilent, Illumina, etc.) to check against
 -type capture_file_type   Type of the capture region file. Allowed file types (default bed): bed, gff and chr-coord (chrXX:coord1-coord2)

** For modifying designs: **
 -modify design_folder     Folder containing the target results of previously calculated capture regions design, which will be modifyed
 -design design_type       Type of design to modify: exon, exon+splicing, exon+UTR, exon+splicing+UTR, design (default)
 -normalize min_size       This option normalizes the intervals so that they have a minimum size of min_size base pairs
 -add intervals_file       Adds the intervals contained in intervals_file to the capture regions in the results_folder
 -subtract intervals_file  Subtracts the intervals contained in intervals_file to the capture regions in the results_folder
 -type design_file_type    Type of the design file. Allowed file types (default bed): bed, gff and chr-coord (chrXX:coord1-coord2)

* Output: *
** In design mode: **
 design/Report.txt              File containing all genes and their bases in exons, splicing sites and UTRs
 design/Log_yyyymmddHHMMSS.txt	File containing all the generated log information
 design/genes/GeneName_all_isoforms.csv, GeneName_exon_coords.txt, GeneName_splicing_coords.txt and GeneName_UTR_coords.txt are the generated files for each individual gene
 design/All_exon_coords.txt, All_exon+splicing_coords.txt, All_exon+splicing+UTR_coords.txt, All_exon+UTR_coords.txt, All_splicing_coords.txt and All_UTR_coords.txt are the files with the coordinates of all the desired transcripts of the given genes in each region. Please note that All_exon+UTR_coords.txt should be the file to use for Agilent's designs

** In check-mode: **
 design/panel_performance/Log.txt
 design/panel_performance/non_covered_bases.txt
 design/panel_performance/Report.txt
 design/panel_performance/genes/GeneName_non_covered_total_bases.txt, GeneName_non_covered_splicing_bases.txt and GeneName_non_covered_UTR_bases.txt

** In modifying designs mode: **
 design/Report_yyyymmddHHMMSS.txt	
 design/Log_yyyymmddHHMMSS.txt	File containing all the generated log information

* Examples of running: *
 \$$scriptname -outdir arbol_mola -prefix BM42_ panel_CSV.txt
 \$$scriptname -target /share/projects/macropanel -check /home/arbol/panel/BaitTiling_bed
 \$$scriptname -modify /share/projects/macropanel -design exon+splicing -add /home/arbol/intervals.gff
 \n";

	exit(1);
}

sub extract_biobase_info{
	# This subroutine extracts the information relating to a variant from biobase info file (previously obtained for a gene)
	
	# Arguments:
	my ($file, $variant, $additional_info) = @_;
	
	# Vars definition:
	my $biobase_info = "";
	
	# Open biobase file for reading
	open BIOBASE, '<', $file or die ("ERROR: Impossible to create $file!! $!");
	
	while (<BIOBASE>){
		my $line = $_;
		my @biobase_fields = split ('\t', $line);
		if ($biobase_fields[0] eq $variant){
			$biobase_info = $biobase_info . $additional_info . $line;
		}
	}
	
	close BIOBASE;
	return $biobase_info;
}

sub gene_coords_parser {
	# This subroutine parses the contents of the output file created by the Biomart query, 
	# checks genes codes & names, duplicates, etc.
	# Finally, it writes the results to a .txt file (tab separated values) that can be uploaded
	# to Agilent's eArray
	# Some other files are created containing:
	# - Splicing zones ($splicing_bases nt at each side of the interval)
	# - UTR zones. In this case, if NO UTR is found for the transcript, $UTR5 and $UTR3 nucleotides are added in the 5' or 3'
	
	# Arguments:
	# $genes_file: csv file containing all the exon coordinates (for all transcripts of the gene)
	# $exon_file: name of the file to output the normalized exon coordinates for desired transcripts
	# $splicing_file: name of the file to output the normalized exon coordinates including splicing sites (30 nt at each side of the exon)
	# $UTR_file: name of the file to output the normalized exon coordinates including splicing sites (30 nt at each side of the exon)
	# $gene_name: name of the gene to consider
	# $transcript: list of the transcripts to consider from gene. All transcripts are considered if empty. This parameter is passed as a reference
	# $sum_exons: variable to return the amount of bases covered by the normalized exon intervals
	# $sum_splicing_exons: variable to return the amount of fases covered by the normalized exon intervals including splicing sites
	# $sum_UTRs: variable to return the amount of fases covered by the UTRs
	my ($genes_file,$LOG,$exon_file,$splicing_file,$strong_splicing_file,$UTR_file,$totals_file,$gene_name,$gene_ensembl_id,$transcript,$sum_exons,$sum_splicing_exons,$sum_strong_splicing_exons,$sum_UTRs,$sum_total,$intervals,$splicing_intervals,$strong_splicing_intervals,$UTR_intervals,$total_intervals,$splicing_bases,$UTR5,$UTR3,$ensembl_UTR5,$ensembl_UTR3,$strong_splicing_bases,$biobase_transcripts,$biobase_UTR5,$biobase_UTR3) = @_;
	
	# variables:
	my $csv = Text::CSV->new();
	my @biomart_row;
	my $row = 0;
	my %exons = ();
	my %totals = ();
	my %splicing = ();
	my %strong_splicing = ();
	my %UTR =();
	my %transcripts = (); 	#Transcripts hash, each $transcript{ENSEMBL_ID} containing:
							#{size}, {UTR5}, {UTR3}, {strand}, {found_UTR5}, {found_UTR3}, {exons}
	my $chromosome = "";
	
	$$gene_ensembl_id = "";
	
	# Open and parse information in $genes_file
	open (CSV, "<", $genes_file) or die("ERROR: Couldn't open file $genes_file!! $!\n");
	
	# Variable for storing the chromosome of the last histocompatibility complex found, so that only a warning from each histocomplex will be displayed
	my $histocomplex_chr = "";
	
	# Count all transcripts present in biomart csv file
	my %all_transcripts;
	while (<CSV>) {
		$row ++;
		if ( $row > 1 ) {
			if ($csv->parse($_)) {
				@biomart_row = $csv->fields();
				$all_transcripts{$biomart_row[1]}=1;
			}
		}
	}
	close CSV;
	
	# Open and parse information in $genes_file
	open (CSV, "<", $genes_file) or die("ERROR: Couldn't open file $genes_file!! $!\n");
	
	$row = 0;
	while (<CSV>) {
		$row ++;
		if ( $row > 1 ) {
			if ($csv->parse($_)) {
				@biomart_row = $csv->fields();
				
				# Check that the gene name and chromosome is the same for all exons
				if ($chromosome eq "" || $chromosome ne $biomart_row[6]) {
					my @chromosomes = (1..22, "X", "Y", "M", "MT");
					foreach (@chromosomes){
						if ( $biomart_row[6] eq $_ ){
							if ($biomart_row[6] eq "MT"){
								$biomart_row[6] = "M";
							}
							$chromosome = $biomart_row[6];
						}
					}
					if ($chromosome eq "" || $chromosome ne $biomart_row[6]){
						if ($histocomplex_chr ne $biomart_row[6]){
							my $line = "WARNING: Found histocompatibility complex in gene $$gene_name. Coordinates from the complex in $biomart_row[6] won't be included in the design!\n";
							print $line;
							print $LOG $line;
							$histocomplex_chr = $biomart_row[6];
						}
					}
				}
				if ($chromosome ne "" && $chromosome eq $biomart_row[6]) {
					if ( length($$gene_name) == 15 && substr($$gene_name,0,4) eq "ENSG" ) { $$gene_name = $biomart_row[7]; }
					if ($biomart_row[7] ne $$gene_name || $biomart_row[6] ne $chromosome ) {
						my $line = "WARNING: Problem in row $row del fichero $genes_file. Gene name ($biomart_row[8]) not expected (different to $$gene_name).\n";
						print $line;
						print $LOG $line;
					} else {
						$$gene_ensembl_id = $biomart_row[0];
	
						# Check all the transcripts and include their exon coordinates in the output file
						my $trans_amount=scalar(@{$transcript});
						my $all_trans = 0;
						if ($trans_amount == 0){
							$trans_amount = scalar(keys(%all_transcripts));
							$all_trans = 1;
						}
						
						# If NO transcript has been defined, then all exon coordinates will be considered.
						# Otherwise, only the exon coordinates from the desired transcripts
						for (my $i=0;$i<$trans_amount;$i++){
							if ( $all_trans == 1 || $biomart_row[1] eq ${$transcript}[$i] ) {
								unless ( exists($transcripts{$biomart_row[1]}) ) {
									$transcripts{$biomart_row[1]}{'size'} = 0;
									if ( $biomart_row[10] == 1 ){
										$transcripts{$biomart_row[1]}{'UTR5'} = min ($biomart_row[8],$biomart_row[9]);
										$transcripts{$biomart_row[1]}{'UTR3'} = max ($biomart_row[8],$biomart_row[9]);
									}
									else {
										$transcripts{$biomart_row[1]}{'UTR5'} = max ($biomart_row[8],$biomart_row[9]);
										$transcripts{$biomart_row[1]}{'UTR3'} = min ($biomart_row[8],$biomart_row[9]);
									}
									$transcripts{$biomart_row[1]}{'strand'} = $biomart_row[10];
									$transcripts{$biomart_row[1]}{'found_UTR5'} = 0;
									$transcripts{$biomart_row[1]}{'found_UTR3'} = 0;   
								}
								
								my @new_interval=($biomart_row[8], $biomart_row[9]);
								
								# Store size and nucleotides of end (position close to coding exon) of the 3'UTR and 5'UTR zones
								$transcripts{$biomart_row[1]}{'size'} = $transcripts{$biomart_row[1]}{'size'} + abs($biomart_row[9] - $biomart_row[8]);
								if ( $transcripts{$biomart_row[1]}{'found_UTR5'} == 0 ) {
									if ( $biomart_row[10] == 1 ) {
										$transcripts{$biomart_row[1]}{'UTR5'} = min ($biomart_row[8],$biomart_row[9],$transcripts{$biomart_row[1]}{'UTR5'});
									}
									else {
										$transcripts{$biomart_row[1]}{'UTR5'} = max ($biomart_row[8],$biomart_row[9],$transcripts{$biomart_row[1]}{'UTR5'});
									}
								}
								if ( $transcripts{$biomart_row[1]}{'found_UTR3'} == 0 ) {
									if ( $biomart_row[10] == 1 ) {
										$transcripts{$biomart_row[1]}{'UTR3'} = max ($biomart_row[8],$biomart_row[9],$transcripts{$biomart_row[1]}{'UTR3'});
									}
									else {
										$transcripts{$biomart_row[1]}{'UTR3'} = min ($biomart_row[8],$biomart_row[9],$transcripts{$biomart_row[1]}{'UTR3'});
									}
								}
								
								# Add intervals to their arrays:
								add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'exons'}}, \@new_interval );
								add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'total_zones'}}, \@new_interval );
								
								# Check if current transcript is a biobase transcript
								my $biobase_transcript = 0;
								if ($#{$biobase_transcripts} > -1){
									foreach (@{$biobase_transcripts}){
										if ($biomart_row[1] eq $_){
											$biobase_transcript = 1;
										}
									}
								}
								
								# Check the 5' UTR from every transcript and include them in the intervals (if corresponding option is set)
								if ( $biomart_row[2] ne "" && $biomart_row[3] ne "" ) {
									my @new_UTR_interval = ($biomart_row[2], $biomart_row[3]); 
									
									# Subtract UTR interval to exons intervals
									subtract_intervals(\@{$transcripts{$biomart_row[1]}{'exons'}},\@new_UTR_interval);
									
									# If ensembl UTRs are to be considered, then its positions are stored in the %transcripts hash, in total_zones and UTRs
									# If not, they get removed from total_zones
									if ( $ensembl_UTR5 == 1 ){
										if ( $biobase_UTR5 == 0 or ( $biobase_UTR5 and $biobase_transcript ) ){
											$transcripts{$biomart_row[1]}{'found_UTR5'} = 1;
											if ( $biomart_row[10] == 1 ) { $transcripts{$biomart_row[1]}{'UTR5'} = min($biomart_row[2],$biomart_row[3]); }
											else { $transcripts{$biomart_row[1]}{'UTR5'} = max($biomart_row[2],$biomart_row[3]); }
											add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'UTRs'}}, \@new_UTR_interval );
											add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'total_zones'}}, \@new_UTR_interval );
										}
									} else {
										subtract_intervals(\@{$transcripts{$biomart_row[1]}{'total_zones'}},\@new_UTR_interval);
									}
								}
								# Check the 3' UTR from every transcript and include them in the intervals
								if ( $biomart_row[4] ne "" && $biomart_row[5] ne "" ) {
									my @new_UTR_interval = ($biomart_row[4], $biomart_row[5]);
									
									# Subtract UTR interval to exons intervals
									subtract_intervals(\@{$transcripts{$biomart_row[1]}{'exons'}},\@new_UTR_interval);
									
									# If ensembl UTRs are to be considered, then its positions are stored in the %transcripts hash, in total_zones and UTRs
									# If not, they get removed from total_zones
									if ( $ensembl_UTR3 == 1 ){
										if ( $biobase_UTR3 == 0 or ( $biobase_UTR3 and $biobase_transcript ) ){
											$transcripts{$biomart_row[1]}{'found_UTR3'} = 1;
											if ( $biomart_row[10] == 1 ) { $transcripts{$biomart_row[1]}{'UTR3'} = max($biomart_row[4],$biomart_row[5]); }
											else { $transcripts{$biomart_row[1]}{'UTR3'} = min($biomart_row[4],$biomart_row[5]); }
											add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'UTRs'}}, \@new_UTR_interval );
											add_and_short_intervals( \@{$transcripts{$biomart_row[1]}{'total_zones'}}, \@new_UTR_interval );
										}
									} else {
										subtract_intervals(\@{$transcripts{$biomart_row[1]}{'total_zones'}},\@new_UTR_interval);
									}
								}
							}
						}
					}
				}
			} else {
				my $err = $csv->error_input;
				print "ERROR: Failed to parse line $row: $err\n";
				print $LOG "ERROR: Failed to parse line $row: $err\n";
			}
		}
	}
	if ($chromosome eq ""){
		print "ERROR: NO coordinates from gene $$gene_name were included in the design!\n";
	}
	
	# The UTR zones will be checked. If not present in none of the chosen transcripts, $UTR5 and $UTR3 nucleotides will be added to the transcripts included in BioBase 
	#(which usually are the interesting ones for mollecular genetists)
	my $one_UTR5_found = 0;
	my $one_UTR3_found = 0;
	foreach ( keys %transcripts ) {
		if ( $transcripts{$_}{'found_UTR5'} != 0 ) {
			$one_UTR5_found++;
		}
		if ( $transcripts{$_}{'found_UTR3'} != 0 ) {
			$one_UTR3_found++;
		}
	}
	if ( $one_UTR5_found == 0 and $UTR5 > 0){
		if ($#{$biobase_transcripts} > -1){
			foreach my $one_transcript (@{$biobase_transcripts}){
				my @new_UTR_interval;
				if (defined($transcripts{$one_transcript})){
					if ( $transcripts{$one_transcript}{'strand'} == 1 ) {
						@new_UTR_interval = ( $transcripts{$one_transcript}{'UTR5'} - $UTR5, $transcripts{$one_transcript}{'UTR5'} );
						add_and_short_intervals ( \@{$transcripts{$one_transcript}{'UTRs'}}, \@new_UTR_interval );
						add_and_short_intervals( \@{$transcripts{$one_transcript}{'total_zones'}}, \@new_UTR_interval );
					} else {
						@new_UTR_interval = ( $transcripts{$one_transcript}{'UTR5'}, $transcripts{$one_transcript}{'UTR5'} + $UTR5 );
						add_and_short_intervals ( \@{$transcripts{$one_transcript}{'UTRs'}}, \@new_UTR_interval );
						add_and_short_intervals( \@{$transcripts{$one_transcript}{'total_zones'}}, \@new_UTR_interval );
					}
					my $line = "INFO: Gene $$gene_name. 5' UTR ( $new_UTR_interval[0],$new_UTR_interval[1] ) will be added to $one_transcript\n";
					print $LOG $line;
					print $line;
				} else {
					my $line = "INFO: HGMD transcript $one_transcript from gene $$gene_name is not included in design\n";
					print $LOG $line;
					print $line;
				}
			}			
		} else{
			# If NO transcripts from the gene are included in Biobase, look for the longest one and add the UTR region
			my $max_size = 0;
			my $biggest_transcript = "";
			foreach ( keys %transcripts ) {
				if ($transcripts{$_}{'size'} > $max_size){
					$biggest_transcript = $_;
					$max_size = $transcripts{$_}{'size'};
				}
			}
			# Only add the UTR if there are existent intervals for the gene
			if ($max_size > 0){
				my @new_UTR_interval;
				if ( $transcripts{$biggest_transcript}{'strand'} == 1 ) {
					@new_UTR_interval = ( $transcripts{$biggest_transcript}{'UTR5'} -$UTR5, $transcripts{$biggest_transcript}{'UTR5'} );
					add_and_short_intervals ( \@{$transcripts{$biggest_transcript}{'UTRs'}}, \@new_UTR_interval );
					add_and_short_intervals( \@{$transcripts{$biggest_transcript}{'total_zones'}}, \@new_UTR_interval );
				} else {
					@new_UTR_interval = ( $transcripts{$biggest_transcript}{'UTR5'}, $transcripts{$biggest_transcript}{'UTR5'} +$UTR5 );
					add_and_short_intervals ( \@{$transcripts{$biggest_transcript}{'UTRs'}}, \@new_UTR_interval );
					add_and_short_intervals( \@{$transcripts{$biggest_transcript}{'total_zones'}}, \@new_UTR_interval );
				}
				my $line = "INFO: No RefSeq transcripts and no 5'-UTRs found for gene $$gene_name. 5' UTR ( $new_UTR_interval[0],$new_UTR_interval[1] ) will be added to $biggest_transcript\n";
				print $LOG $line;
				print $line;
			}
		}
	}
	if ( $one_UTR3_found == 0 and $UTR3 > 0){
		if ($#{$biobase_transcripts} > -1){
			foreach my $one_transcript (@{$biobase_transcripts}){
				my @new_UTR_interval;
				if (defined($transcripts{$one_transcript})){
					if ( $transcripts{$one_transcript}{'strand'} == 1 ) {
						@new_UTR_interval = ( $transcripts{$one_transcript}{'UTR3'}, $transcripts{$one_transcript}{'UTR3'} + $UTR3 );
						add_and_short_intervals ( \@{$transcripts{$one_transcript}{'UTRs'}}, \@new_UTR_interval );
						add_and_short_intervals( \@{$transcripts{$one_transcript}{'total_zones'}}, \@new_UTR_interval );
					} else {
						@new_UTR_interval = ( $transcripts{$one_transcript}{'UTR3'} - $UTR3, $transcripts{$one_transcript}{'UTR3'} );
						add_and_short_intervals ( \@{$transcripts{$one_transcript}{'UTRs'}}, \@new_UTR_interval );
						add_and_short_intervals( \@{$transcripts{$one_transcript}{'total_zones'}}, \@new_UTR_interval );
					}
					my $line = "INFO: Gene $$gene_name. 3' UTR ( $new_UTR_interval[0],$new_UTR_interval[1] ) will be added to $one_transcript\n";
					print $LOG $line;
					print $line;
				} else {
					my $line = "INFO: HGMD transcript $one_transcript from gene $$gene_name is not included in design\n";
					print $LOG $line;
					print $line;
				}
			}
		} else{
			# If NO transcripts from the gene are included in Biobase, look for the longest one and add the UTR region
			my $max_size = 0;
			my $biggest_transcript = "";
			foreach ( keys %transcripts ) {
				if ($transcripts{$_}{'size'} > $max_size){
					$biggest_transcript = $_;
					$max_size = $transcripts{$_}{'size'};
				}
			}
			# Only add the UTR if there are existent intervals for the gene
			if ($max_size > 0){
				my @new_UTR_interval;
				if ( $transcripts{$biggest_transcript}{'strand'} == 1 ) {
					@new_UTR_interval = ( $transcripts{$biggest_transcript}{'UTR3'}, $transcripts{$biggest_transcript}{'UTR3'} + $UTR3 );
					add_and_short_intervals ( \@{$transcripts{$biggest_transcript}{'UTRs'}}, \@new_UTR_interval );
					add_and_short_intervals( \@{$transcripts{$biggest_transcript}{'total_zones'}}, \@new_UTR_interval );
				} else {
					@new_UTR_interval = ( $transcripts{$biggest_transcript}{'UTR3'} - $UTR3, $transcripts{$biggest_transcript}{'UTR3'} );
					add_and_short_intervals ( \@{$transcripts{$biggest_transcript}{'UTRs'}}, \@new_UTR_interval );
					add_and_short_intervals( \@{$transcripts{$biggest_transcript}{'total_zones'}}, \@new_UTR_interval );
				}
				my $line = "INFO: No RefSeq transcripts and no 3'-UTRs found for gene $$gene_name. 3' UTR ( $new_UTR_interval[0],$new_UTR_interval[1] ) will be added to $biggest_transcript\n";
				print $LOG $line;
				print $line;
			}
		}
	}
	
	# Creation and update of the intervals hashes, which contain all the intervals of all the chromosomes involved in the design
	for (keys %transcripts){
		my $one_transcript = $_;
		add_intervals_to_hash ( $intervals, \@{$transcripts{$one_transcript}{'exons'}}, $chromosome );
		add_intervals_to_hash ( \%exons, \@{$transcripts{$one_transcript}{'exons'}}, $chromosome );
		
		# Obtain splicing intervals and add them to the global and local splicing intervals hashes
		my $splicing_temp = \@{ obtain_splicing_from_exons(\@{$transcripts{$one_transcript}{'exons'}}, $splicing_bases) };
		add_intervals_to_hash ( $splicing_intervals, $splicing_temp, $chromosome );
		add_intervals_to_hash ( \%splicing, $splicing_temp, $chromosome );
		$splicing_temp = \@{ obtain_splicing_from_exons(\@{$transcripts{$one_transcript}{'UTRs'}}, $splicing_bases) };
		add_intervals_to_hash ( $splicing_intervals, $splicing_temp, $chromosome );
		add_intervals_to_hash ( \%splicing, $splicing_temp, $chromosome );

		# Obtain strong splicing intervals and add them to the global and local strong splicing intervals hashes
		$splicing_temp = \@{ obtain_splicing_from_exons(\@{$transcripts{$one_transcript}{'exons'}}, $strong_splicing_bases) };
		add_intervals_to_hash ( $strong_splicing_intervals, $splicing_temp, $chromosome );
		add_intervals_to_hash ( \%strong_splicing, $splicing_temp, $chromosome );
		$splicing_temp = \@{ obtain_splicing_from_exons(\@{$transcripts{$one_transcript}{'UTRs'}}, $strong_splicing_bases) };
		add_intervals_to_hash ( $strong_splicing_intervals, $splicing_temp, $chromosome );
		add_intervals_to_hash ( \%strong_splicing, $splicing_temp, $chromosome );
		
		# Add UTRs to the global and local UTRs intervals hashes		
		add_intervals_to_hash ( $UTR_intervals, \@{$transcripts{$one_transcript}{'UTRs'}}, $chromosome );
		add_intervals_to_hash ( \%UTR, \@{$transcripts{$one_transcript}{'UTRs'}}, $chromosome );
		
		# Add the splicing intervals to the total target intervals:
		add_intervals_to_hash ( $total_intervals, \@{$transcripts{$one_transcript}{'total_zones'}}, $chromosome );
		add_intervals_to_hash ( \%totals, \@{$transcripts{$one_transcript}{'total_zones'}}, $chromosome );
	}
	
	# Close file
	close CSV;
	
	# Chromosome name "MT" is renamed to "M"
	if ( $chromosome eq "MT" ) { $chromosome = "M"; }
	
	# Count and print intervals
	$$sum_exons = print_and_count_intervals(\@{$exons{$chromosome}},$chromosome,$exon_file,0);
	$$sum_total = print_and_count_intervals(\@{$totals{$chromosome}},$chromosome,$totals_file,0);
	$$sum_splicing_exons = print_and_count_intervals(\@{$splicing{$chromosome}},$chromosome,$splicing_file,0);
	$$sum_strong_splicing_exons = print_and_count_intervals(\@{$strong_splicing{$chromosome}},$chromosome,$strong_splicing_file,0);
	$$sum_UTRs = print_and_count_intervals(\@{$UTR{$chromosome}},$chromosome,$UTR_file,0);
	
	# Display error if no gene with the given name was found!
	if ( $row < 2 ) {
		print "ERROR: Problem with gene $$gene_name!! No information found!!\n";
		print $LOG "ERROR: Problem with gene $$gene_name!! No information found!!\n";
		$chromosome = "ERROR";
	}
	
	return $chromosome;
}

sub refseq_transcripts_parser {
	# This subroutine parses the contents of the output refseq-mRNAs file created by the Biomart query. 
	# It checks the existent refseq transcripts for the given genes, and it retrieves the biobase-variants information for these transcripts belonging to a gene
	# (if present in HGMD -biobase-)
	# Finally, it creates an intervals hash with following information:
	# intervals_hash: {chr}{transcript}{variant}[[coordinate_begin,coordinate_final]]
	# It also stores the variants information into a gene related file
	# It returns the number of Biobase variants found for the gene
	
	# Arguments:
	my ($refseq_file,$biobase_vars_file,$LOG,$biobase_files,$vars_file,$gene,$biobase_transcripts) = @_;
	
	# variables:
	my @chroms = (1..22,"X","Y");
	my %vars_intervals = ();

	# Open output file for writing:
	open OUTPUT, '>', $biobase_vars_file or die ("ERROR: Impossible to create $biobase_vars_file!! $!");
	
	# Take each refseq entry and check if its present in the biobase files
	open (REFSEQ, "<", $refseq_file) or die("ERROR: Couldn't open file $refseq_file!! $!\n");
	my $found_vars = 0;
	my $row_refseq = 0;
	while (<REFSEQ>){
		$row_refseq ++;
		if ($row_refseq > 1){
			my $line = $_;
			chomp $line;
			my @cols = split(',', $line);
			
			# Store the biobase transcripts (refseq IDs):
			if (defined($cols[1])){
				my $refseq = $cols[1];
				push(@{$biobase_transcripts}, $cols[0]);
			
				# Look for the refseq entry in the different chromosome-biobase files
				my $i=0;
				while($i<$#chroms){
					# Open and parse information in biobase files
					my $biobase_file = $biobase_files . $chroms[$i] . ".txt";
					open (BIOBASE, "<", $biobase_file) or die("ERROR: Couldn't open file $biobase_file!! $!\n");
					while (<BIOBASE>) {
						my @biobase_fields = split ('\t', $_);
						my @variant = split ('\.', $biobase_fields[9]);
						if ( $variant[0] eq $refseq ) {
							my @values = split(':', $biobase_fields[4]);
							my @tmp = split('chr',$values[0]);
							my @coords = split('-',$values[1]);
	
							# Detect problems in biobase
							if (! defined($coords[1]) or $coords[1] eq "" or $coords[1] eq "-"){
								$coords[1] = $coords[0];
							}
							
							# Store the coordinates of the variants in a intervals hash:
							$found_vars ++;
							my @new_interval = [($coords[0],$coords[1])];
							add_intervals_to_hash(\%vars_intervals,\@new_interval,$chroms[$i],$biobase_fields[1]);
													
							# Output to gene_biobase file:
							print OUTPUT "$biobase_fields[1]\t$biobase_fields[7]\t$biobase_fields[8]\t$biobase_fields[9]\thttp://www.ncbi.nlm.nih.gov/pubmed?term=$biobase_fields[11]\n";
						}
					}
					close BIOBASE;
					
					# Increment chromosome counter:
					$i++;
				}
			}
		}
	}
	if ($row_refseq < 2){
		print $LOG "WARNING: No RefSeq equivalent found for gene ${gene}!!\n";
		print "WARNING: No RefSeq equivalent found for gene ${gene}!!\n";
	}
	if ($found_vars == 0){
		print $LOG "WARNING: No Biobase entries found for gene ${gene}!!\n";
		print "WARNING: No Biobase entries found for gene ${gene}!!\n";
	}
	close REFSEQ;
	close OUTPUT;
	
	my $vars_bases = create_intervals_file ($vars_file, "eArray", $LOG, \%vars_intervals);
	return $row_refseq - 1;
}

sub biomart_query{
	# Function for retrieving information from biomart queries to ensembl
	# $output_type: is either CSV or TSV
	# $query_type: is either gene_info or refseq_info
	
	# Arguments:
	my ( $gene_name, $gene, $query_type, $output_file, $output_type, $confFile, $LOG, $attempt, $attempts) = @_;
	
	# An example script demonstrating the use of BioMart API.
	# This perl API representation is only available for configuration versions >=  0.5 
	
	#
	# NB: change action to 'clean' if you wish to start a fresh configuration  
	# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
	#
	
	my $action='cached';
	my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
	my $registry = $initializer->getRegistry;
	
	my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
		
		$query->setDataset("hsapiens_gene_ensembl");
		if ( length($gene_name)==15 && substr($gene_name,0,4) eq "ENSG" ) { $query->addFilter("ensembl_gene_id", @{$gene}); }
		else { $query->addFilter("hgnc_symbol", @{$gene}); }
		if ($query_type eq "gene_info"){
			#$query->addFilter("biotype", ["protein_coding"]);
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
			$query->addAttribute("strand");
		} elsif ($query_type eq "refseq_info"){
			#$query->addFilter("biotype", ["protein_coding"]);
			$query->addAttribute("ensembl_transcript_id");
			$query->addAttribute("refseq_mrna");
		}
	
	$query->formatter($output_type);
	
	my $query_runner = BioMart::QueryRunner->new();
	
	# Redirect STDOUT to output_file:
	{
		# Open the output file for the results (it will capture STDOUT)
		local *STDOUT;
		open(STDOUT, '>', $output_file) || die("ERROR: Can't redirect STDOUT to $output_file: $!");
		
		############################## GET RESULTS ##########################
		# to obtain unique rows only
		$query_runner->uniqueRowsOnly(1);
		
		$query_runner->execute($query);
		$query_runner->printHeader();
		$query_runner->printResults();
		$query_runner->printFooter();
		#####################################################################

		# close the file for outputting the results:
		close STDOUT;
	}
	
	# Check output file. If empty, retry the query until max_try!!
	open(FILE, '<', $output_file) || die("ERROR: Can't check $output_file: $!");
	my @file = <FILE>;
	close FILE;
	if ($#{@file} == 0){
		if ($attempt < $attempts){
			$attempt ++;
			my $line = "INFO: Possible problem in connection to Ensembl-biomart. Retrying connection! Attempt $attempt\n";
			print $LOG $line;
			print $line;
			biomart_query( $gene_name, $gene, $query_type, $output_file, $output_type, $confFile, $LOG, $attempt, $attempts);
		}
	}
}

# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
our $fooGlobalVar = $ENV{"PATH"};
my $strong_splicing_bases = 5; #Variable temporarily hardcoded... sorry!!

# ------------------------
# -- default parameters --
# ------------------------
my $prefix = "";
my $outdir = cwd() . "/panel_info";
my $genes_file = "";
my $transcript = "";
my $splicing_bases = 20;
my $minsize = 0;
my $UTR5 = 158;
my $UTR3 = 502;
my $ensembl_UTR5 = 1;
my $ensembl_UTR3 = 1;
my $biobase_UTR5 = 0;
my $biobase_UTR3 = 0;
my $run_flag = 0;
my $target_folder = "";
my $panel_design = "";
my $file_type = "bed";
my $biomart_conf = "/share/apps/src/biomart-perl/conf/ensembl.xml";
my $file_intervals_to_add = "";
my $modify_folder = "";
my $design_type = "";
my $normalize_minsize = "";
my $add_file = "";
my $subtract_file = "";
my $biobase_files = "/share/references/biobase/bioBaseSortUniques_chr";
my $attempts = 3;

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v1.4
Script for obtaining gene coordinates from the Ensembl-biomart databases, to check its results against designs (Agilent, Illumina...) and to modify previously created designs
:: Copyright 2012 by Arbol :: 
******************************************************************************************************************************************************************************\n";

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
        elsif ($opt eq "-splicing") {
			$splicing_bases = shift @ARGV;
        }
        elsif ($opt eq "-minsize") {
			$minsize = shift @ARGV;
        }
        elsif ($opt eq "-UTR5") {
			$UTR5 = shift @ARGV;
        }
        elsif ($opt eq "-UTR3") {
			$UTR3 = shift @ARGV;
        }
        elsif ($opt eq "-no_ensembl_UTR5") {
			$ensembl_UTR5 = 0;
        }
		elsif ($opt eq "-no_ensembl_UTR3") {
			$ensembl_UTR3 = 0;
        }
        elsif ($opt eq "-only_biobase_UTR5") {
			$biobase_UTR5 = 1;
        }
        elsif ($opt eq "-only_biobase_UTR3") {
			$biobase_UTR3 = 1;
        }
		elsif ($opt eq "-target") {
			$target_folder = shift @ARGV;
			$run_flag = 1;
        }
        elsif ($opt eq "-check") {
			$panel_design = shift @ARGV;
			$run_flag = 1;
        }
		elsif ($opt eq "-biomart") {
			$biomart_conf = shift @ARGV;
        }
		elsif ($opt eq "-type") {
			$file_type = shift @ARGV;
        }
        elsif ($opt eq "-intervals") {
			$file_intervals_to_add = shift @ARGV;
        }
		elsif ($opt eq "-modify") {
			$modify_folder = shift @ARGV;
			$run_flag = 2;
        }
        elsif ($opt eq "-design") {
			$design_type = shift @ARGV;
        }
        elsif ($opt eq "-normalize") {
			$normalize_minsize = shift @ARGV;
        }
        elsif ($opt eq "-add") {
			$add_file = shift @ARGV;
        }
		elsif ($opt eq "-subtract") {
			$subtract_file = shift @ARGV;
        }
        elsif ($opt eq "-attempts") {
			$attempts = shift @ARGV;
        }
		else {
			print "ERROR: Bad option: $opt\n";
			print "Call $basename$ext with no parameters for help.\n";
			exit(1);
		}
	}
	else {
		# Mandatory parameters definitions
		$genes_file = $opt;
	}
}
if ($genes_file eq "" && $run_flag == 0) {
	print "ERROR: Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ($run_flag == 1 && ($target_folder eq "" || $panel_design eq "")) {
	print "ERROR: Please be sure to include both -target and -check parameters for check mode in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if ($run_flag ==2 && ($design_type ne "design" && $design_type ne "exon" && $design_type ne "exon+splicing" && $design_type ne "exon+UTR" && $design_type ne "exon+splicing+UTR")){
	print "ERROR: Please be sure to include an allowed design type (design,exon,exon+splicing,exon+UTR,exon+splicing+UTR) in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);	
}
if ($file_type ne "bed" && $file_type ne "gff" && $file_type ne "chr-coord" && $file_type ne "tsv"){
	print "ERROR: Unsupported design coordinates file type!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

#####################################################################################################
# Biomart data retrieval mode - design mode -
#####################################################################################################
if ($run_flag == 0){
	# Create output folders:
	unless(-d $outdir){
		mkdir $outdir or die ("ERROR: Impossible to create folder $outdir!!");
	}
	unless(-d $outdir . "/genes"){
		mkdir $outdir . "/genes" or die ("ERROR: Impossible to create folder $outdir/genes!!");
	}
	
	# Create and open log file:
	my $yyyymmdd = strftime("%Y%m%d%H%M%S", localtime);
	my $log_file = $outdir . "/". $prefix . "Log_" . $yyyymmdd . ".txt";
	open my $LOG, '>', $log_file or die ("ERROR: Impossible to create $log_file!! $!");
	
	my $now = localtime;
	print "INFO: Script started at $now\n";
	print $LOG "INFO: Script started at $now\n";
	
	# Vars definition
	my $csv = Text::CSV->new(); # Variable for dealing with the CSV file
	my @file_fields; # Array containing the fields of the input file
	my @genes_array; # Array containing rows with gene name followed by desired transcripts
	my @genes_ensembl; # Array containing the ensembl IDs of the genes in the panel design
	my %intervals = (); # Hash of exon intervals
	my %splicing_intervals = (); # Hash of splicing sites intervals
	my %strong_splicing_intervals = (); # Hash of splicing sites intervals
	my %UTR_intervals = (); # Hash of UTR intervals
	my %total_intervals = (); # Hash of total intervals (exons, splicing sites, UTRs)
	
	# Parse the input CSV files with genes names and transcripts
	open (CSV, "<", $genes_file) or die $!;
	my $row = 0;
	while (<CSV>) {
		if ($csv->parse($_)) {
			@file_fields = $csv->fields();
			print "Gene: $file_fields[0] Transcripts: @file_fields[1..$#file_fields]\n";
			print $LOG "Gene: $file_fields[0] Transcripts: @file_fields[1..$#file_fields]\n";
			push (@{$genes_array[$row]}, @file_fields);
			$row ++;
		} else {
			my $err = $csv->error_input;
			print "ERROR: Failed to parse line: $err";
		}
	}
	close CSV;
	
	# Create, open and print header of stats file:
	my $report_file=$outdir . "/". "Report.txt";
	open STATS, '>', $report_file or die ("ERROR: Impossible to create $report_file!! $!");
	my $line = "Gene symbol (HGNC)\tENSEMBL Gene ID\tChromosome\tExon bases\tStrong splicing sites bases\tSplicing sites bases\tUTRs bases\tTotal bases\n";
	print STATS $line;
	
	# Obtain coordinates and parameters for each gene
	for ( my $i=0; $i<=$#genes_array; $i++ ) {
		
		# The whitespaces from the beginning and end of the gene names should be removed:
		my $gene_name = $genes_array[$i][0];
		my $gene_ensembl_id;
		
		# Biomart needs an array as input:
		$gene_name =~ s/^\s+//;
		$gene_name =~ s/\s+$//;
		my @gene = [$gene_name];
		
		my $output_file = $outdir . "/genes/". $gene_name . "_all_isoforms.csv";
		my $refseq_file = $outdir . "/genes/". $gene_name . "_refseq_mRNA.csv";
		print "\nRetrieving results for gene $gene_name to file $output_file\n";
		print $LOG "\nRetrieving results for gene $gene_name to file $output_file\n";
		
		# Queries to  Ensembl via biomart (these queries are separated since otherwise biomart returns an error due to the amount of external queries):
		biomart_query($gene_name,\@gene,"gene_info",$output_file,"CSV",$biomart_conf,$LOG,0,$attempts);
		biomart_query($gene_name,\@gene,"refseq_info",$refseq_file,"CSV",$biomart_conf,$LOG,0,$attempts);
		
		# Define output files and variables:
		my $cols = @{ $genes_array[$i] };
		my @transcripts = @{ $genes_array[$i] } [ 1..($cols-1) ];
		my $exons_file = $outdir . "/genes/" . $gene_name . "_exon_coords.txt";
		my $splicing_file = $outdir . "/genes/" . $gene_name . "_splicing_coords.txt";
		my $strong_splicing_file = $outdir . "/genes/" . $gene_name . "_strong_splicing_coords.txt";
		my $UTR_file = $outdir . "/genes/" . $gene_name . "_UTR_coords.txt";
		my $totals_file = $outdir . "/genes/" . $gene_name . "_totals_coords.txt";
		my $biobase_file = $outdir . "/genes/" . $gene_name . "_biobase_vars.txt";
		my $vars_file = $outdir . "/genes/" . $gene_name . "_vars_coords.txt";
		my $sum_desired_exons;
		my $sum_desired_splicing_exons;
		my $sum_strong_splicing_exons;
		my $sum_UTRs;
		my $sum_total;
		my @biobase_transcripts;
		
		# Parse the refseq files and check the biobase files to retrieve interesting variants
		my $vars_amount = refseq_transcripts_parser($refseq_file,$biobase_file,$LOG,$biobase_files,$vars_file,$gene_name,\@biobase_transcripts);
		
		# Parse the generated gene files, compute results and store in intervals hashes
		# Please note that in the call to gene_coords_parser, some parameters are passed as references!!
		my $chromosome = gene_coords_parser( $output_file, $LOG, $exons_file, $splicing_file, $strong_splicing_file, $UTR_file, $totals_file, \$gene_name, \$gene_ensembl_id, \@transcripts, \$sum_desired_exons, \$sum_desired_splicing_exons, \$sum_strong_splicing_exons, \$sum_UTRs, \$sum_total, \%intervals, \%splicing_intervals, \%strong_splicing_intervals, \%UTR_intervals, \%total_intervals, $splicing_bases, $UTR5, $UTR3, $ensembl_UTR5, $ensembl_UTR3, $strong_splicing_bases, \@biobase_transcripts, $biobase_UTR5, $biobase_UTR3 );		
		
		# Check if the gene was already in the list of genes in the panel design (it's possible that some synonim genes have been included in the design)
		# If the gene hadn't been yet included in the list, then it will be included
		my $gene_found = 0;
		foreach (@genes_ensembl) {
			if ($_ eq $gene_ensembl_id) {$gene_found = 1;}
		}
		if ($gene_found == 0 && $gene_ensembl_id ne "") {
			push(@genes_ensembl,$gene_ensembl_id);
			# Output to report file of gene	0
			$line = $gene_name . "\t" . $gene_ensembl_id . "\t" . $chromosome . "\t" . $sum_desired_exons . "\t" . $sum_strong_splicing_exons . "\t" . $sum_desired_splicing_exons . "\t" . $sum_UTRs . "\t" . $sum_total . "\n";
			print STATS $line;
		} elsif ($gene_ensembl_id ne "") {
			$line = "WARNING: Gene $gene_name ($gene_ensembl_id) was already included in the design!!\n";
			print $line;
			print $LOG $line;
		} else {
			$line = "ERROR: Problem with gene $gene_name ($gene_ensembl_id)!! No information found!!\n";
			print $line;
			print $LOG $line;
		}
		
	}
		
	# Create all the different intervals hashes:
	my %exon_splicing = %{ add_intervals_hashes ( \%intervals, \%splicing_intervals) };
	my %exon_UTRs = %{ add_intervals_hashes ( \%intervals, \%UTR_intervals) };
	my %exon_splicing_UTRs = %{ add_intervals_hashes ( \%exon_splicing, \%UTR_intervals) };
	my %design;
	my %all_design;
	my %intervals_to_add;
	
	my $added_intervals_bases = 0;
	if ( $file_intervals_to_add ne "" ){
		import_intervals_file($file_intervals_to_add, "tsv", \%intervals_to_add, $LOG);
		%design = %{ add_intervals_hashes ( \%exon_UTRs, \%intervals_to_add) };
		%all_design = %{ add_intervals_hashes ( \%exon_splicing_UTRs, \%intervals_to_add ) };
		my $file = $outdir . "/" . "All_added_intervals_coords";
		$added_intervals_bases = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%intervals_to_add);
		$added_intervals_bases = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%intervals_to_add);
		$added_intervals_bases = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%intervals_to_add);		
		%all_design = %{ add_intervals_hashes ( \%exon_splicing_UTRs, \%intervals_to_add ) };		
	} else {
		%design = %exon_UTRs;
		%all_design = %exon_splicing_UTRs;
	}
	
	# Create output files and count bases:
	my $file = $outdir . "/" . "All_exon_coords";
	my $all_bases = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%intervals);
	$all_bases = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%intervals);
	$all_bases = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%intervals);
	
	$file = $outdir . "/" . "All_strong_splicing_coords";
	my $all_strong_splicing = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%strong_splicing_intervals);
	$all_strong_splicing = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%strong_splicing_intervals);
	$all_strong_splicing = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%strong_splicing_intervals);	

	$file = $outdir . "/" . "All_splicing_coords";
	my $all_splicing = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%splicing_intervals);
	$all_splicing = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%splicing_intervals);
	$all_splicing = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%splicing_intervals);
	
	$file = $outdir . "/" . "All_UTR_coords";
	my $all_UTRs = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%UTR_intervals);
	$all_UTRs = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%UTR_intervals);
	$all_UTRs = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%UTR_intervals);
	
	$file = $outdir . "/" . "All_exon+splicing_coords";
	my $all_bases_splicing = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%exon_splicing);
	$all_bases_splicing = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%exon_splicing);
	$all_bases_splicing = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%exon_splicing);
	
	$file = $outdir . "/" . "All_exon+UTR_coords";
	my $all_bases_UTRs = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%exon_UTRs);
	$all_bases_UTRs = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%exon_UTRs);
	$all_bases_UTRs = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%exon_UTRs);
	
	$file = $outdir . "/" . "All_exon+splicing+UTR_coords";
	my $all_bases_splicing_UTRs = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%exon_splicing_UTRs);
	$all_bases_splicing_UTRs = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%exon_splicing_UTRs);
	$all_bases_splicing_UTRs = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%exon_splicing_UTRs);
	
	$file = $outdir . "/" . "Design_coords";
	normalize_intervals_hash (\%design, $LOG, $minsize);
	my $bases_design = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%design);
	$bases_design = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%design);
	$bases_design = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%design);

	$file = $outdir . "/" . "All_design_coords";
	my $all_bases_design = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%all_design);
	$all_bases_design = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%all_design);
	$all_bases_design = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%all_design);
	
	$file = $outdir . "/" . "All_added_intervals";
	$added_intervals_bases = create_intervals_file ($file . "_eArray.txt", "eArray", $LOG, \%intervals_to_add);
	$added_intervals_bases = create_intervals_file ($file . "_TruSeq.txt", "TruSeq", $LOG, \%intervals_to_add);
	$added_intervals_bases = create_intervals_file ($file . "_HaloPlex.txt", "HaloPlex", $LOG, \%intervals_to_add);
	
	# Finish LOG file:
	$line = "------------------------------------------------------------------- \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in exons from target regions: $all_bases bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in strong splicing sites from target regions: $all_strong_splicing bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in splicing sites from target regions: $all_splicing bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in UTRs from target regions: $all_UTRs bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in exons & splicing sites from target regions: $all_bases_splicing bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in exons & splicing sites & UTRs from target regions: $all_bases_splicing_UTRs bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases in exons & UTRs from target regions: $all_bases_UTRs bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases from additional custom intervals: $added_intervals_bases bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Total number of bases from target regions (exons, splicing sites, UTRs and intervals to be added): $all_bases_design bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	$line = "Number of bases of the design to submit (exons, UTRs and intervals to be added): $bases_design bp \n";
	print $line;
	print $LOG $line;
	print STATS $line;
	close STATS;
	
	# Output finishing time in the log file
	$now = localtime;
	print "INFO: Script finished successfully at $now\n";
	print $LOG "INFO: Script finished successfully at $now\n";
	close $LOG;
}

#####################################################################################################
# Check mode!!!!!!!!
##################################################################################################### 
elsif ($run_flag == 1) {
	# Variables definitions:
	my %earray_intervals = ();
	my %design_intervals = ();
	my %exon_intervals = ();
	my %strong_splicing_intervals = ();
	my %splicing_intervals = ();
	my %UTR_intervals = ();
	my %vars_intervals = ();
	my %added_intervals = ();
	my $design_file = $target_folder . "/All_design_coords_eArray.txt";
	my $exon_file = $target_folder . "/All_exon_coords_eArray.txt";
	my $strong_splicing_file = $target_folder . "/All_strong_splicing_coords_eArray.txt";
	my $splicing_file = $target_folder . "/All_splicing_coords_eArray.txt";
	my $UTR_file = $target_folder . "/All_UTR_coords_eArray.txt";
	my $added_intervals_file = $target_folder."/All_added_intervals_eArray.txt";
	
	# Create output folders:
	my $yyyymmdd = strftime("%Y%m%d%H%M%S", localtime);
	my $output_folder = $target_folder . "/" . $prefix . "panel_performance_" . $yyyymmdd;
	unless(-d $output_folder){
		mkdir $output_folder or die ("ERROR: Impossible to create folder $output_folder!!"); 
	}
	my $output_genes_folder = $output_folder . "/genes";
	unless(-d $output_genes_folder){
		mkdir $output_genes_folder or die ("ERROR: Impossible to create folder $output_genes_folder!!"); 
	}
	
	# Create and open log file:
	my $log_file=$output_folder . "/Log.txt";
	open my $LOG, '>', $log_file or die ("ERROR: Impossible to create $log_file!! $!");
	
	my $now = localtime;
	print "INFO: Script started at $now\n";
	print $LOG "INFO: Script started at $now\n";
	
	# Create, open and print header of stats files:
	my $report_file = $output_folder . "/Report.txt";
	open STATS, '>', $report_file or die ("ERROR: Impossible to create $report_file!! $!");
	my $line = "Gene symbol (HGNC)\tENSEMBL Gene ID\tChromosome\tExon bases\tNon-covered exon bases\tPercentage of non-covered exon bases\tStrong splicing sites bases\tNon-covered strong splicing sites bases\tPercentage of non-covered strong splicing site bases\tSplicing sites bases\tNon-covered splicing sites bases\tPercentage of non-covered splicing site bases\tUTRs bases\tNon-covered UTR bases\tPercentage of non-covered UTR bases\tTotal bases\tNon-covered total bases\tPercentage of non-covered total bases\tNumber of Biobase variants\tNon-covered Biobase variants\tPercentage of non-covered Biobase variants\n";
	print STATS $line;
	my $report_vars_file = $output_folder . "/Biobase_variants_not_covered.txt";
	open VARS, '>', $report_vars_file or die ("ERROR: Impossible to create $report_vars_file!! $!");
	$line = "Gene symbol (HGNC)\tENSEMBL Gene ID\tChromosome\tBegin position\tFinal position\tVariant name\tDisease info\tBases close to variant\tRefSeq identification\tPubMed link\n";
	print VARS $line;
	
	# Data imports:
	my $panel_design_bases = import_intervals_file($panel_design, $file_type, \%earray_intervals, $LOG);
	$line = "Number of bases in the imported panel design (SureSelect, HaloPlex, TruSeq, ...): $panel_design_bases\n";
	print $LOG $line;
	print $line;
	
	my $design_bases = import_intervals_file($design_file, "chr-coord", \%design_intervals, $LOG);
	$line = "Number of bases in target regions: $design_bases\n";
	print $LOG $line;
	print $line;
	
	$design_bases = import_intervals_file($exon_file, "chr-coord", \%exon_intervals, $LOG);
	$line = "Number of bases in exons from target regions: $design_bases\n";
	print $LOG $line;
	print $line;

	$design_bases = import_intervals_file($strong_splicing_file, "chr-coord", \%strong_splicing_intervals, $LOG);
	$line = "Number of bases in strong splicing sites from target regions: $design_bases\n";
	print $LOG $line;
	print $line;
		
	$design_bases = import_intervals_file($splicing_file, "chr-coord", \%splicing_intervals, $LOG);
	$line = "Number of bases in splicing sites from target regions: $design_bases\n";
	print $LOG $line;
	print $line;
	
	$design_bases = import_intervals_file($UTR_file, "chr-coord", \%UTR_intervals, $LOG);
	$line = "Number of bases in UTR from target regions: $design_bases\n";
	print $LOG $line;
	print $line;
	
	$design_bases = import_intervals_file($added_intervals_file, "chr-coord", \%added_intervals, $LOG);
	$line = "Number of bases in added intervals: $design_bases\n";
	print $LOG $line;
	print $line;
	
	# Subtract the earray_intervals to the design_intervals:
	my %non_covered_intervals = %{ subtract_intervals_hashes(\%design_intervals,\%earray_intervals) };
	my %non_covered_exon_intervals = %{ subtract_intervals_hashes(\%exon_intervals,\%earray_intervals) };
	my %non_covered_strong_splicing_intervals = %{ subtract_intervals_hashes(\%strong_splicing_intervals,\%earray_intervals) };
	my %non_covered_splicing_intervals = %{ subtract_intervals_hashes(\%splicing_intervals,\%earray_intervals) };
	my %non_covered_UTR_intervals = %{ subtract_intervals_hashes(\%UTR_intervals,\%earray_intervals) };
	my %non_covered_added_intervals = %{ subtract_intervals_hashes(\%added_intervals,\%earray_intervals) };
	
	# Create output files and count bases:
	# All non-covered bases:
	my $file = $output_folder . "/Non_covered_bases.txt";
	my $non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_intervals);
	print $LOG "Number of non-covered bases from the design: $non_covered_bases\n";
	print "Number of non-covered bases from the design: $non_covered_bases\n";
	
	# All non-covered exon bases:
	$file = $output_folder . "/Exons_non_covered.txt";
	$non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_exon_intervals);
	print $LOG "Number of non-covered EXON bases from the design: $non_covered_bases\n";
	print "Number of non-covered EXON bases from the design: $non_covered_bases\n";
	
	# All non-covered strong splicing sites:
	$file = $output_folder . "/Strong_splicing_non_covered.txt";
	$non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_strong_splicing_intervals);
	print $LOG "Number of non-covered STRONG SPLICING bases from the design: $non_covered_bases\n";
	print "Number of non-covered STRONG SPLICING bases from the design: $non_covered_bases\n";

	# All non-covered splicing sites:
	$file = $output_folder . "/Splicing_non_covered.txt";
	$non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_splicing_intervals);
	print $LOG "Number of non-covered SPLICING bases from the design: $non_covered_bases\n";
	print "Number of non-covered SPLICING bases from the design: $non_covered_bases\n";

	# All non-covered UTR bases:
	$file = $output_folder . "/UTR_non_covered.txt";
	$non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_UTR_intervals);
	print $LOG "Number of non-covered UTR bases from the design: $non_covered_bases\n";
	print "Number of non-covered UTR bases from the design: $non_covered_bases\n";
	
	# All non-covered UTR bases:
	$file = $output_folder . "/Added_intervals_non_covered.txt";
	$non_covered_bases = create_intervals_file ($file, "eArray", $LOG, \%non_covered_added_intervals);
	print $LOG "Number of non-covered added intervals bases: $non_covered_bases\n";
	print "Number of non-covered added intervals bases: $non_covered_bases\n";
	
	# Check coverage on each individual gene, using info from design report file:
	$file = $target_folder . "/Report.txt";
	open DESIGN, '<', $file or die ("ERROR: Impossible to open $file!! $!");
	while (<DESIGN>){
		my @design = split('\t',$_);
		if (defined($design[1]) && substr($design[1],0,4) eq "ENSG"){
			# Vars definitions:
			my $gene;
			my %total_intervals = ();
			my %strong_splicing_intervals = ();			
			my %splicing_intervals = ();
			my %UTR_intervals = ();
			my %exon_intervals = ();
			my %vars_intervals = ();
			my ($percentage_total,$percentage_strong_splicing,$percentage_splicing,$percentage_UTR,$percentage_exon,$percentage_vars) = (0,0,0,0,0,0);
			
			opendir(DIR,$target_folder . "/genes") or die "ERROR: Couldn't open directory, $!";
			my $ensembl_id = 0;
			foreach (grep(/$design[1]/,readdir(DIR))) {$ensembl_id=1;}	
			if ($ensembl_id) {$gene = $design[1];}
			else {$gene = $design[0];}
			
			# Import and count intervals files:
			my $total_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_totals_coords.txt", "chr-coord", \%total_intervals, $LOG);
			my $strong_splicing_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_strong_splicing_coords.txt", "chr-coord", \%strong_splicing_intervals, $LOG);
			my $splicing_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_splicing_coords.txt", "chr-coord", \%splicing_intervals, $LOG);
			my $UTR_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_UTR_coords.txt", "chr-coord", \%UTR_intervals, $LOG);
			my $exon_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_exon_coords.txt", "chr-coord", \%exon_intervals, $LOG);
			my $vars_bases = import_intervals_file($target_folder . "/genes/" . $gene . "_vars_coords.txt", "chr-coord", \%vars_intervals, $LOG);
			my $amount_vars = 0;
			foreach (keys %vars_intervals){
				$amount_vars = $amount_vars + scalar(keys %{$vars_intervals{$_}});
			}
			
			# Subtract design intervals and calculate non_covered_bases
			my %non_covered_total_intervals = %{ subtract_intervals_hashes(\%total_intervals,\%earray_intervals) };
			my %non_covered_strong_splicing_intervals = %{ subtract_intervals_hashes(\%strong_splicing_intervals,\%earray_intervals) };
			my %non_covered_splicing_intervals = %{ subtract_intervals_hashes(\%splicing_intervals,\%earray_intervals) };
			my %non_covered_UTR_intervals = %{ subtract_intervals_hashes(\%UTR_intervals,\%earray_intervals) };
			my %non_covered_exon_intervals = %{ subtract_intervals_hashes(\%exon_intervals,\%earray_intervals) };
			my %non_covered_vars_intervals = %{ subtract_intervals_hashes(\%vars_intervals,\%earray_intervals) };
			
			# Calculate amount of non covered bases and output non covered intervals to file
			my $non_covered_total_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_total_bases.txt", "eArray", $LOG, \%non_covered_total_intervals);
			if ($total_bases != 0){$percentage_total = $non_covered_total_bases / $total_bases;}
			my $non_covered_strong_splicing_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_strong_splicing_bases.txt", "eArray", $LOG, \%non_covered_strong_splicing_intervals);
			if ($splicing_bases != 0){$percentage_strong_splicing = $non_covered_strong_splicing_bases / $strong_splicing_bases;}
			my $non_covered_splicing_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_splicing_bases.txt", "eArray", $LOG, \%non_covered_splicing_intervals);
			if ($splicing_bases != 0){$percentage_splicing = $non_covered_splicing_bases / $splicing_bases;}
			my $non_covered_UTR_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_UTR_bases.txt", "eArray", $LOG, \%non_covered_UTR_intervals);
			if ($UTR_bases != 0){$percentage_UTR = $non_covered_UTR_bases / $UTR_bases;}
			my $non_covered_exon_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_exon_bases.txt", "eArray", $LOG, \%non_covered_exon_intervals);
			if ($exon_bases != 0){$percentage_exon = $non_covered_exon_bases / $exon_bases;}
			my $non_covered_vars_bases = create_intervals_file ($output_genes_folder . "/" . $gene . "_non_covered_vars_bases.txt", "HaloPlex", $LOG, \%non_covered_vars_intervals);
	
			# Output to BioBase and report files
			my $amount_non_covered_vars = 0;
			foreach ( keys %non_covered_vars_intervals ){
				my $chr = $_;
				$amount_non_covered_vars = $amount_non_covered_vars + scalar(keys %{$non_covered_vars_intervals{$chr}});
				foreach ( keys %{$non_covered_vars_intervals{$chr}} ){
					my $var = $_;
					foreach ( @{$non_covered_vars_intervals{$chr}{$var}}) {
						$line = extract_biobase_info($target_folder . "/genes/" . $gene . "_biobase_vars.txt", $var, $design[0] . "\t" . $design[1] . "\t" . $chr . "\t" . ${$_}[0] . "\t" . ${$_}[1] . "\t");
					}
					print VARS $line;
				}
			}
			if ($amount_vars != 0){$percentage_vars = $amount_non_covered_vars / $amount_vars;}
			print STATS "$design[0]\t$design[1]\t$design[2]\t$exon_bases\t$non_covered_exon_bases\t$percentage_exon\t$strong_splicing_bases\t$non_covered_strong_splicing_bases\t$percentage_strong_splicing\t$splicing_bases\t$non_covered_splicing_bases\t$percentage_splicing\t$UTR_bases\t$non_covered_UTR_bases\t$percentage_UTR\t$total_bases\t$non_covered_total_bases\t$percentage_total\t$amount_vars\t$amount_non_covered_vars\t$percentage_vars\n";
			$line = "Number of non-covered bases from gene ${gene}: $non_covered_total_bases\n";
			print $LOG $line;
			print $line;			
		}
	}
	close DESIGN;
	close STATS;
	close VARS;
	
	# Output finishing time in the log file
	$now = localtime;
	print "INFO: Script finished successfully at $now\n";
	print $LOG "INFO: Script finished successfully at $now\n";
	close $LOG;
}


#####################################################################################################
# Modify mode!!!!!!!!
##################################################################################################### 
else {
	# Variables definitions:
	my %design_intervals = ();
	my $design_file;
	
	# Create and open log file:
	my $yyyymmdd = strftime("%Y%m%d%H%M%S", localtime);
	my $log_file = $modify_folder . "/". $prefix . "Log_" . $yyyymmdd . ".txt";
	open my $LOG, '>', $log_file or die ("ERROR: Impossible to create $log_file!! $!");
	
	# Open for appending the stats file:
	my $report_file = $modify_folder . "/". $prefix . "Report_" . $yyyymmdd . ".txt";
	open STATS, '>>', $report_file or die ("ERROR: Impossible to open for appending $report_file!! $!");
	
	# Display start time:
	my $now = localtime;
	print "INFO: Script started at $now\n";
	print $LOG "INFO: Script started at $now\n";
	
	# Design to be modified:
	if ($design_type eq "design"){
		$design_file = $modify_folder . "/" . $prefix . "Design_coords_";	
	} else {
		$design_file = $modify_folder . "/" . $prefix . "All_" . $design_type . "_coords_";
	}
	
	# Data imports:
	import_intervals_file($design_file . "eArray.txt", "chr-coord", \%design_intervals, $LOG);
	
	# Backup original files:
	my @filetypes = ("eArray","TruSeq","HaloPlex");
	foreach (@filetypes){
		my $filetobecopied = $design_file . $_ . ".txt";
		my $newfile = $design_file . $_ . "_" . $yyyymmdd . ".txt";
		copy($filetobecopied, $newfile) or die "ERROR: File $filetobecopied cannot be copied!!\n";	
	}
	
	# Normalize intervals:
	if ($normalize_minsize ne ""){
		my $line = "INFO: Normalizing intervals files with names starting with $design_file\n";
		print STATS $line;
		print $LOG $line;
		print $line;
		normalize_intervals_hash (\%design_intervals, $LOG, $normalize_minsize);
	}
	
	# Add intervals:
	if ($add_file ne ""){
		my $line = "INFO: Adding file $add_file to intervals files with names starting with $design_file\n";
		print STATS $line;
		print $LOG $line;
		print $line;
		my %intervals_to_add = ();
		import_intervals_file($add_file, $file_type, \%intervals_to_add, $LOG);
		%design_intervals = %{ add_intervals_hashes (\%design_intervals,\%intervals_to_add)};
	}
	
	# Subtract intervals:
	if ($subtract_file ne ""){
		my $line = "INFO: Subtracting file $add_file to intervals files with names starting with $design_file\n";
		print STATS $line;
		print $LOG $line;
		print $line;
		my %intervals_to_subtract = ();
		import_intervals_file($subtract_file, $file_type, \%intervals_to_subtract, $LOG);
		%design_intervals = %{ subtract_intervals_hashes (\%design_intervals,\%intervals_to_subtract)};		
	}
	
	# Save changes to files:
	foreach (@filetypes){
		my $total_bases = create_intervals_file ($design_file . $_ . ".txt", $_, $LOG, \%design_intervals);
		my $line = "The number of base-pairs in file $design_file$_.txt has being changed to: $total_bases\n";
		print STATS $line;
		print $LOG $line;
		print $line;
	}
	
	# Output finishing time in the log file
	$now = localtime;
	print "INFO: Script finished successfully at $now\n";
	print $LOG "INFO: Script finished successfully at $now\n";
	close $LOG;
	close STATS;
}