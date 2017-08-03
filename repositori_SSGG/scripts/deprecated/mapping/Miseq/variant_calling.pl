#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to identify variants
# @Author: JM Rosa
# @Contributors: Arbol
###############################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use file_handle;
use json_handle;
use NGS::Illumina::essay;
use NGS::Illumina::var_calling;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $metadata;
my $hold_jid = 1;
my $cwd = `pwd`;
chomp $cwd;
my $output_id = cwd() . "/output_file";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"m=s"	=> \$metadata,
				"o=s"	=> \$output_id,
				"hold_jid=s"	=> \$hold_jid,
);

print "\n=================================================================
$basename: Script to call and filter variants v1.0\n";

if (not defined $metadata) {
	die "
			Options:
		-m Metadata File
		-o Ouput name (optional)
		-hold_jid (Job Identifier to wait for initiation, optional)
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	print "Running $basename with parameters m=$metadata hold_jid=$hold_jid and o=$output_id\n";
	
	my $essay = essay->get_essay($metadata, $cwd);
	my $var_id = var_calling->variant_calling ($essay, $cwd, $hold_jid);
	print "JOB_ID\t$var_id\n";
}


exit;