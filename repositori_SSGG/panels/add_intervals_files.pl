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
use genomics::intervals qw(add_intervals_arrays subtract_intervals_arrays show_intervals_array add_and_short_intervals subtract_intervals count_intervals print_and_count_intervals show_intervals_hash add_intervals_to_hash add_intervals_hashes subtract_intervals_from_hash subtract_intervals_hashes count_intervals_hash print_and_count_intervals_hash normalize_intervals_hash obtain_splicing_from_exons);
use genomics::transcripts qw(new_transcript show_transcript);
use utils::bioinfo_files_utilities qw(import_intervals_file create_intervals_file);

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for adding two intervals files (currently bed, gff, gff3, tsv or chr-coord files supported)

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 file1 file2		Files containing the intervals to be added

* Options: *
 -o output_file     Name of the output file. Default output is sum.bed
 -t1 file1_type     Type of file1 (bed, gff, gff3, tsv or chr-coord). Default is bed
 -t2 file1_type     Type of file2 (bed, gff, gff3, tsv or chr-coord). Default is bed
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
my $output = "sum.bed";
my $type1 = "bed";
my $type2 = "bed";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.1
Script for adding two intervals files (currently bed, gff, gff3, tsv or chr-coord files supported)
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
		if ($opt eq "-o") {
			$output = shift @ARGV;
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
		$file2 = shift @ARGV;
	}
}

if ($file1 eq "" or $file2 eq "") {
	print "ERROR: Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}

#################################################################################################3
# ---------------
# -- main body -- 
# ---------------
################################################################################################3

# Variables definitions:
my %file1_intervals = ();
my %file2_intervals = ();
	
# Display start time:
my $now = localtime;
print "INFO: Script started at $now\n";
	
# Data imports:
print "INFO: Importing $file1...\n";
import_intervals_file($file1, $type1, \%file1_intervals, *STDOUT);
#show_intervals_hash(\%file1_intervals);
print "INFO: Importing $file2...\n";
import_intervals_file($file2, $type2, \%file2_intervals, *STDOUT);
#show_intervals_hash(\%file2_intervals);

# Add intervals:
my $line = "INFO: Adding intervals file $file1 to $file2\n";
print $line;
my %result_intervals = %{ add_intervals_hashes (\%file1_intervals,\%file2_intervals)};
	
# Save changes to file
my $total_bases = create_intervals_file ($output, "bed", *STDOUT, \%result_intervals);
$line = "The number of base-pairs in sum file $output is $total_bases\n";
print $line;
	
# Output finishing time in the log file
$now = localtime;
print "INFO: Script finished successfully at $now\n";
#******************************************************************* 
