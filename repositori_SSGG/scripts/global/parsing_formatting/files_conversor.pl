#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
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
use Text::CSV;
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];
use POSIX qw( strftime );
use File::Copy;
use utils::bioinfo_files_utilities qw(import_intervals_file create_intervals_file file_conversion);

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for converting intervals files from a given format to a different one (currently bed, gff, gff3, tsv or chr-coord files supported)

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 input_file            File to be converted

* Options: *
 -o output_file        Name of the output file. Default output is sum.bed
 -t1 input_file_type   Type of input file (bed, gff, gff3, tsv or chr-coord). Default is chr-coord
 -t2 output_file_type  Type of output file (bed, gff, gff3, tsv or chr-coord). Default is bed
  \n";

	exit(1);
}


# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------

# ------------------------
# -- default parameters --
# ------------------------
my $prefix = "";
my $file1 = "";
my $file2 = "";
my $type1 = "chr-coord";
my $type2 = "bed";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.1
Script for converting intervals files from a given format to a different one (currently bed, gff, gff3, tsv or chr-coord files supported)
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
	my $opt = shift @ARGV;
	if ($opt =~ /^\-/) {
		if ($opt eq "-o") {
			$file2 = shift @ARGV;
        } elsif ($opt eq "-t1") {
			$type1 = shift @ARGV;
        } elsif ($opt eq "-t2") {
			$type2 = shift @ARGV;
        } else {
			print "ERROR: Bad option: $opt\n";
			print "Call $basename$ext with no parameters for help.\n";
			exit(1);
		}
	}
	else {
		# Mandatory parameters definitions
		$file1 = $opt;
	}
}

if ($file1 eq "") {
	print "ERROR: Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

# Display start time:
my $now = localtime;
print "INFO: Script started at $now\n";

if ($file2 eq ""){
	(my $output_basename, my $output_path, my $input_ext) = fileparse($file1, qr/\.[^.]*/);
	$file2 = $output_path . $output_basename . "." . $type2;
}
	
# File conversion:
file_conversion ($file1, $type1, $file2, $type2, *STDOUT);
	
# Output finishing time in the log file
$now = localtime;
print "INFO: Script finished successfully at $now\n";
#******************************************************************* 
