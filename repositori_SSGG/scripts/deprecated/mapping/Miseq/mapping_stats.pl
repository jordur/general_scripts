#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to map Illumina data with BWA
# @Author: Arbol
# @Contributors: JM Rosa
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
use NGS::mapping_stats;
use NGS::essay;

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
my $offset = 1;
my $coverages = "1,10,20,50,100,200,400,600";
my $slots = 8;
my $queue = "shudra";
my $id = "";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"m=s"	=> \$metadata,
				"coverages=s"	=> \$coverages,
				"offset"	=> \$offset,
				"o=s"	=> \$output_id,
				"slots=i"	=> $slots,
				"q=s"	=> $queue,
				"hold_jid=s"	=> $id,
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
	print Dumper($essay);
	exit;
	my $stats_id = mapping_stats->stats_essay($essay, $cwd);
	print "JOB_ID\t$stats_id\n";
}

exit;