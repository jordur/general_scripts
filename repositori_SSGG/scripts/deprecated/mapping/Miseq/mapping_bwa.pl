#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to map Illumina data with BWA
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
use NGS::Illumina::mapping_bwa;
use NGS::Illumina::essay;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $metadata;
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
);

print "\n=================================================================
$basename: Script to map Illumina data with BWA v1.0\n";

if (not defined $metadata) {
	die "
			Options:
		-m Metadata File
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "Process finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	print "Running $basename with parameters m=$metadata and o=$output_id\n";
	
	my $essay = essay->get_essay($metadata, $cwd);
	my $aln_id = mapping_bwa->mapping_essay ($essay, $cwd);
	print "JOB_ID\t$aln_id\n";
}


exit;