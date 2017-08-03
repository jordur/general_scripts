#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for launching pipelines to analyze
# @Author: Arbol
# @Contributors: Arbol
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use File::Basename;
use Cwd;
use Cwd qw( abs_path );
use strict;
use warnings;
use POSIX qw( strftime );
use JSON -support_by_pp;
use Data::Dumper;
use Getopt::Long;

# ---------------
# -- functions --
# ---------------

sub usage ();

# -----------------------
# -- global parameters --
# -----------------------
#my $fooGlobalVar = $ENV{"PATH"};

# ------------------------
# -- default parameters --
# ------------------------
my $conf_file = "/share/apps/etc/pipelines/config.json";

&main();

sub main (){
	# Get input parameters
	GetOptions (	
					"c=s"	=> \$conf_file,
	);
	
	# Display name and version, and check input arguments
	(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);
	print "$basename$ext v0.1
	Script for launching analysis pipelines (PIPElines To Analyze->PipeTA). It loads configurations, paths and folders, so that it assures the use of the right information. · Copyright 2012 by Arbol 
	***************************************************************************************************************************************************\n";
	if (not defined $conf_file) { usage; }
	
	# Load default configuration from XML file:
	my $json;
	{
		local $/; #Enable 'slurp' mode
		open my $fh, "<", $conf_file;
		$json = <$fh>;
		close $fh;
	}
	my $data = decode_json($json);
	print Dumper($data); #uncomment for showing contents of configuration stored in $data
}

sub usage () {
# Displays the script usage
	# Parameters:
	my ( $scriptname ) = @_;

	print "
Usage:
  $scriptname <options>

* Input (mandatory parameters): *

* Options: *
-c Configuration_file    File containing all the configuration parameters that must be loaded (default /share/apps/etc/pipelines/config.json)
* Output: *

* Examples of running: *
 \$$scriptname 
 \n";
 
	exit(1);
}

exit;