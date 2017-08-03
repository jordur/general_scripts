#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script for transforming intervals files from different formats
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
use Getopt::Long;
use Data::Dumper;
use Text::CSV;
use bioinfo_files_utilities qw(import_intervals_file create_intervals_file);
use intervals qw(subtract_intervals_hashes add_intervals_hashes);

# ---------------------------
# -- functions definitions --
# ---------------------------


# -----------------------
# -- global parameters --
# -----------------------
my $input;
my $output;
my $inputformat = "chr-coord";
my $outputformat = "bed";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"i=s"	=> \$input,
				"o=s"	=> \$output,
				"ti=s"	=> \$inputformat,
				"to=s"	=> \$outputformat,
);

print "\n=================================================================
$basename: is a script for transforming intervals files from different formats (currently bed, gff, gff3, tsv or chr-coord files supported) v0.1\n";

if (not defined $input) {
	die "
* Input (mandatory parameters): *
 -i input_file           Intervals file that is going to be transformed

* Options: *
 -o output_file          Name of the output file. Default output is input.ext
 -ti input file format   Format of input file (bed, gff, gff3, tsv or chr-coord). Default is bed
 -to output file format  Format of output file (bed, gff, gff3, tsv or chr-coord). Default is bed
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main(){
	my $path = `pwd`;
	chomp $path;
	my %input_intervals;
	
	(my $basename, my $dir, my $ext) = fileparse($input, qr/\.[^.]*/);
	if (not defined($output)){
		$output = $basename.".".$outputformat;
	}
	
	# Data imports:
	print "INFO: Importing $input...\n";
	import_intervals_file($input, $inputformat, \%input_intervals, *STDOUT);
	
	# Save changes to file
	my $total_bases = create_intervals_file ($output, $outputformat, *STDOUT, \%input_intervals);
	my $line = "The number of base-pairs in input file $input is $total_bases\n";
	print $line;
}

exit;