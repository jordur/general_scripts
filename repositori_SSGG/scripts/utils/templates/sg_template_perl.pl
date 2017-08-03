#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $input;
my $output = cwd() . "/output_file.tsv";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"f=s"	=> \$input,
				"o=s"	=> \$output,
);

print "\n=================================================================
$basename: Template script for developing Perl scripts v1.1\n";

if (not defined $input) {
	die "
			Options:
		-f chr-coord file (chrXX:coord1-coord2 files)
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	my $path = `pwd`;
	chomp $path;
	
	print "Running $basename with parameters f=$input, o=$output\n";
	
}

exit;