#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for joining two bam/sam files in a single one
# @Author: Arbol
# @Contributors: JMRosa
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use File::Basename;
use Cwd;

# ---------------
# -- functions --
# ---------------

sub usage {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
$scriptname is a script for joining two bam/sam files in a single one.
Basically, the script adds '_F3' at the end of the F3 ID_RD, and changes the ID_SAMPLE to a common name
(alternatives: samtools view -h name.bam -Shb)
(picard tools for adding infos)

Usage:
  $scriptname <options> -paired input1 -F3 input2

* Input (mandatory parameters): *
 -paired input1     csfasta file
 -F3 input2			qual file

* Options: *
 -prefix sampleprefix_      Adds the given string as prefix for the ouput results
 -outdir output_folder      Creates the output files in the given output folder, creating it if needed

* Examples of running: *
 \$$scriptname -prefix ohio_ -outdir arbol_mola -paired PE.bam -F3 F3.bam
 \n";

	exit(1);
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
my $outdir = cwd();
my $paired = "";
my $F3 = "";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v1.0
Script for joining two bam/sam files in a single one · Copyright 2012 by Arbol 
******************************************************************************\n";

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
		elsif ($opt eq "-paired") {
			$paired = shift @ARGV;
        }
		elsif ($opt eq "-F3") {
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
	}
}
if ($paired eq "" or $F3 eq "") {
	print "Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}


# ---------------
# -- main body -- 
# ---------------

my $now = localtime;
print "Script started at $now\n";

# Check if the files are bam or sam and add the needed options
(my $paired_name, my $paired_dir, my $paired_ext) = fileparse($paired, qr/\.[^.]*/);
(my $F3_name, my $F3_dir, my $F3_ext) = fileparse($F3, qr/\.[^.]*/);
if ($paired_ext eq ".sam") {$options_paired = $options_paired . " -S";}
if ($F3_ext eq ".sam") {$options_F3 = $options_F3 . " -S";}

# Parse the pair-end file:
my $command = `samtools view $options_paired`;

$now = localtime;
print "Script finished at $now\n";

#******************************************************************* 
