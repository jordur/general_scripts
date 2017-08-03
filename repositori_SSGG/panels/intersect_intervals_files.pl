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
use Text::CSV;
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];
use POSIX qw( strftime );
use File::Copy;
use bioinfo_files_utilities qw(import_intervals_file create_intervals_file);
use intervals qw(subtract_intervals_hashes add_intervals_hashes);

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for intersecting two intervals files (currently bed, gff, gff3, tsv or chr-coord files supported)

Usage:
  $scriptname <options> input

* Input (mandatory parameters): *
 file1 file2		Files containing the intervals to be intersected

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
my $output = "intersection.bed";
my $type1 = "bed";
my $type2 = "bed";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v0.1
Script for adding two intervals files (currently bed, gff, gff3, tsv or chr-coord files supported)
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

# Intersect intervals:
my $line = "INFO: Intersecting intervals files $file1 and $file2\n";
print $line;

my %temp_intervals = %{ subtract_intervals_hashes (\%file1_intervals,\%file2_intervals)};
my %result_intervals = %{ subtract_intervals_hashes (\%file1_intervals,\%temp_intervals)};

# Save changes to file
my $total_bases = create_intervals_file ($output, "bed", *STDOUT, \%result_intervals);
$line = "The number of base-pairs in intersection file $output is $total_bases\n";
print $line;
	
# Output finishing time in the log file
$now = localtime;
print "INFO: Script finished successfully at $now\n";