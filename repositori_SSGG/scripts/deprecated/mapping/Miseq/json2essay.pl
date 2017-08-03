#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to create essay for pipelines
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
use NGS::Illumina::json2essay;
use file_handle;
use json_handle;

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
$basename: Script to create essays for pipelines v1.0\n";

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
	
	my %data = &get_from_json($metadata);
	
	my @essays = json2essay->get_essays(\%data, $cwd);
	
	my $output = &get_out_file_handle($output_id, 'tsv');
	print $output "#Essay\tjson_file\taligner\n";
	
	foreach my $essay (@essays){
		my $json = &create_json ($essay);
		print $output "$essay->{name}\t$json\t$essay->{aligner}\n";
	}
}

exit;