#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Arbol
# @Contributors:
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
$scriptname is a template script for developing Perl scripts.

Usage:
  $scriptname <options> input1 input2

* Input (mandatory parameters): *
 input1     csfasta file
 input2     qual file

* Options: *
 -prefix sampleprefix_      Adds the given string as prefix for the ouput results
 -outdir output_folder      Creates the output files in the given output folder, creating it if needed

* Examples of running: *
 \$$scriptname -prefix ohio_ file1.csfasta file2_QV.qual
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
my $csfasta = "";
my $qual = "";

# -----------------------------
# -- Name, version and usage --
# -----------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
print "$basename$ext v1.0
Template script for developing of Perl scripts · Copyright 2012 by Arbol 
************************************************************************\n";

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
		$csfasta = $opt;
		$qual = shift @ARGV;
	}
}
if ($csfasta eq "" or $qual eq "") {
	print "Not enough parameters in call to $basename$ext!!\n";
	print "Call $basename$ext with no parameters for help.\n";
	exit(1);
}


# ---------------
# -- main body -- 
# ---------------

$now = localtime;
print "Script started at $now\n";

# Add here the script code | | | | | |
#                          V V V V V V
print "Hi there! Command launched: $0 prefix=$prefix outdir=$outdir csfasta=$csfasta qual=$qual \n";

$now = localtime;
print "Script finished at $now\n";

#******************************************************************* 
