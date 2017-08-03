#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for filtering out Rfam miRNAs and repeated miRNAs in closer regions
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
use Text::CSV;
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];
use POSIX qw( strftime );
use File::Copy;
use lib '/share/apps/scripts/global/SSGG_libraries/genomics';
use intervals qw(add_intervals_arrays subtract_intervals_arrays show_intervals_array add_and_short_intervals subtract_intervals count_intervals print_and_count_intervals show_intervals_hash add_intervals_to_hash add_intervals_hashes subtract_intervals_from_hash subtract_intervals_hashes count_intervals_hash print_and_count_intervals_hash normalize_intervals_hash overlapping_intervals_hashes filtering_repeats);
use transcripts qw(new_transcript show_transcript);
use bioinfo_files_utilities qw(import_intervals_file create_intervals_file import_miRDeep2_file);
#use Data::Dumper;

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for filtering out Rfam miRNAs and repeated miRNAs in closer regions from a list of novel miRNAs 

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 input      File containing the miRNAs obtained (default format is miRDeep2 csv file)

* Options: *
 -itype file_type       Type of input file. Allowed file types (default miRDeep2): miRDeep2, bed, gff, gff3 and chr (chrXX:coord1-coord2)
 -otype file_type       Type of input file. Allowed file types (default bed): bed, gff, gff3 and chr (chrXX:coord1-coord2)
 -prefix sampleprefix_  Adds the given string as prefix for the output results
 -Rfam path_to_Rfam     Path to the folder containing the Rfam miRNAs (default /share/references/lifescope_referenceData/lifetech/hg19/Rfam)
 -species_files list    List of the Rfam files to consider, those belonging to the studied specie sorted according to chromosome. Default are Homo Sapiens files: CM000663.1.gff3,CM000664.1.gff3,CM000665.1.gff3,CM000666.1.gff3,CM000667.1.gff3,CM000668.1.gff3,CM000669.1.gff3,CM000670.1.gff3,CM000671.1.gff3,CM000672.1.gff3,CM000673.1.gff3,CM000674.1.gff3,CM000675.1.gff3,CM000676.1.gff3,CM000677.1.gff3,CM000678.1.gff3,CM000679.1.gff3,CM000680.1.gff3,CM000681.1.gff3,CM000682.1.gff3,CM000683.1.gff3,CM000684.1.gff3,CM000685.1.gff3,CM000686.1.gff3
 -Rfam_overlap #bases   Maximum allowed percentage of overlapping bases between new miRNAs and Rfam miRNAs (default 50%)
 -repeat_overlap #bases Maximum allowed number of overlapping bases between new miRNAs in repeated zones (default 1bp)

* Output: *
 prefix_input_filtered  File containing the filtered miRNAs

* Examples of running: *
 \$$scriptname -outdir arbol_mola -prefix BM42_ panel_CSV.txt
 \n";
 
	exit(1);
}

sub extract_biobase_info{
	# This subroutine extracts the information relating to a variant from biobase info file (previously obtained for a gene)
	
	# Arguments:
	my ($file, $variant) = @_;
	
	# Vars definition:
	my $biobase_info = "";
	
	# Open biobase file for reading
	open BIOBASE, '<', $file or die ("ERROR: Impossible to create $file!! $!");
	
	while (<BIOBASE>){
		my $line = $_;
		my @biobase_fields = split ('\t', $line);
		if ($biobase_fields[0] eq $variant){
			$biobase_info = $biobase_info . $line;
		}
	}
	
	close BIOBASE;
	return $biobase_info;
}


# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
our $fooGlobalVar = $ENV{"PATH"};

# ------------------------
# -- default parameters --
# ------------------------
my $miRNAs_file = "";
my $prefix = "";
my $outdir = cwd();
my $itype = "miRDeep2";
my $otype = "bed";
my $Rfam = "/share/references/lifescope_referenceData/lifetech/hg19/Rfam";
my @species_files = ( "CM000663.1.gff3", "CM000664.1.gff3", "CM000665.1.gff3", "CM000666.1.gff3", "CM000667.1.gff3", "CM000668.1.gff3", "CM000669.1.gff3", "CM000670.1.gff3", "CM000671.1.gff3", "CM000672.1.gff3", "CM000673.1.gff3", "CM000674.1.gff3", "CM000675.1.gff3", "CM000676.1.gff3", "CM000677.1.gff3", "CM000678.1.gff3", "CM000679.1.gff3", "CM000680.1.gff3", "CM000681.1.gff3", "CM000682.1.gff3", "CM000683.1.gff3", "CM000684.1.gff3", "CM000685.1.gff3", "CM000686.1.gff3");
my $Rfam_overlap = 50;
my $repeat_overlap = 1; #bases Maximum allowed number of overlapping bases between new miRNAs in repeated zones (default 1bp)

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.1
Script for filtering out Rfam miRNAs and repeated miRNAs in closer regions from a list of novel miRNAs
· Copyright 2012 by Arbol · 
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
		if ($opt eq "-prefix") {
			$prefix = shift @ARGV;
        }
		elsif ($opt eq "-itype") {
			$itype = shift @ARGV;
        }
        elsif ($opt eq "-otype") {
			$otype = shift @ARGV;
        }
        elsif ($opt eq "-Rfam") {
			$Rfam = shift @ARGV;
        }
		elsif ($opt eq "-species_files") {
			my $species_list = shift @ARGV;
			@species_files = split(',', $species_list);
        }
        elsif ($opt eq "-Rfam_overlap") {
			$Rfam_overlap = shift @ARGV;
        }
        elsif ($opt eq "-repeat_overlap") {
			$repeat_overlap = shift @ARGV;
        }
		else {
			print "ERROR: Bad option: $opt\n";
			print "Call $basename$ext with no parameters for help.\n";
			exit(1);
		}
	}
	else {
		# Mandatory parameters definitions
		$miRNAs_file = $opt;
	}
}
if ($miRNAs_file eq ""){
	print "ERROR: Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}
if (($itype ne "bed" and $itype ne "gff" and $itype ne "chr" and $itype ne "miRDeep2") or ($otype ne "bed" and $otype ne "gff" and $otype ne "chr")){
	print "ERROR: Not allowed file type!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

# Output file:
my $output_Rfam = $prefix . $miRNAs_file . "_Rfam_filtered." . $otype;
my $output_file = $prefix . $miRNAs_file . "_filtered." . $otype;

# Create and open log file:
my $yyyymmdd = strftime("%Y%m%d%H%M%S", localtime);
my $log_file = $outdir . "/". $prefix . "Log_" . $yyyymmdd . ".txt";
open my $LOG, '>', $log_file or die ("ERROR: Impossible to create $log_file!! $!");

my $now = localtime;
print "INFO: Script started at $now\n";
print $LOG "INFO: Script started at $now\n";

# Vars definition
my %miRNAs_hash = ();

# Import input file
import_miRDeep2_file($miRNAs_file, \%miRNAs_hash, $LOG);
#show_intervals_hash(\%miRNAs_hash);

# Filter out Rfam miRNAs
# Check all the Rfams:
print "Starting filtering out Rfam miRNAs...\n";
foreach (@species_files){
	print "Parsing contents of Rfam $_\n";
	my %Rfam_hash = ();
	
	# Import Rfam miRNAs
	import_intervals_file($Rfam . "/" . $_, "gff3", \%Rfam_hash, $LOG);
	
	# Check each miRNA interval. If it intersects a given Rfam, dismish it
	overlapping_intervals_hashes(\%miRNAs_hash,\%Rfam_hash,$Rfam_overlap,$repeat_overlap,"ge",$LOG);
}
# Output to file
create_intervals_file($output_Rfam,"HaloPlex",$LOG,\%miRNAs_hash);

# Filter out repeated miRNAs
# Check each miRNA with all the miRNAs in the file. Dismish those that are intersecting more than $repeat_overlap bp
print "Starting filtering out repeated miRNAs...\n";
filtering_repeats(\%miRNAs_hash,$repeat_overlap,$LOG);

# Output to file
create_intervals_file($output_file,"HaloPlex",$LOG,\%miRNAs_hash);