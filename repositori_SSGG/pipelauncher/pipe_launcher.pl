#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to launch Miseq Pipeline with BWA
# @Author: Arbol & Rosa
# @Contributors: 
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
use NGS::Illumina::mapping_bwa;
use NGS::Illumina::mapping_novoalign;
use NGS::SOLiD::mapping_bioscope;
use NGS::primary_stats;
use NGS::mapping_stats;
use NGS::var_calling;
use NGS::annotation;
use NGS::project;
use file_handle;
use json_handle;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $metadata;
my $queue = "dalits.q";
my $priority = 0;
my $cwd = `pwd`;
chomp $cwd;

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"m=s"	=> \$metadata,
				"q=s"	=> \$queue,
				"p=i"	=> \$priority,
);

print "\n=================================================================
$basename: PipeTA: Pipeline to Analyze, a launcher for in silico analysis from metadata information... v1.3\n";

if (not defined $metadata) {
	die "
			Options:
		-m Metadata File
		-q queue to run analysis on
		-p priority of generated jobs
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "Process finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){
	print "Running $basename with parameters of $metadata, in queue $queue with priority $priority\n";
	
	# Create project object from json file
	my $project = project->get_project($metadata,$queue,$priority);
	
	# Check all essays and launch their modules
	foreach (keys(%{$project->{essays}})){
		my $essay = $project->{essays}{$_};
		print "\nLaunching jobs for $essay->{name} ($essay->{seq_platform} platform)\n";
		# Launch analysis modules

		# Do not care of getting jobs_ids that the new module must wait for, since this task
		# will be carried out in submit.pm (with info from state.json)
		foreach my $module (@{$essay->{flow_list}}){
			my $state = substr($essay->get_parameter("modules",$module,"state"),0,8);
			my $job_id = ""; #$essay->GetModuleHoldJobs($module);
			if ($state ne "finished"){
				if ($module eq "mapping"){
					my $aligner = $essay->get_parameter("modules","mapping","parameters","aligner");
					if ($essay->is_solid() eq 'true'){
						if ($aligner eq "bioscope"){
							mapping_bioscope->mapping($essay, $job_id);
						}
					} else {
						if ($aligner eq "bwa"){
							mapping_bwa->mapping($essay, $job_id);
						} elsif ($aligner eq "novoalign"){
							mapping_novoalign->mapping($essay, $job_id);
						}
					}
				} elsif ($module eq "primary_stats"){
					primary_stats->stats($essay,$job_id);
				} elsif ($module eq "mapping_stats"){
					mapping_stats->stats($essay,$job_id);
				} elsif ($module eq "variant_calling"){
					var_calling->variant_calling($essay, $job_id);
				} elsif ($module eq "annotation"){
					annotation->annotate($essay, $job_id);
				} elsif ($module eq "validation"){
					validation->validate($essay,$job_id);
				}
				# Update essay state file
				create_json($essay,$essay->{path} . "/state");
			}
		}
		
		# Update essay state file
		create_json($essay,$essay->{path} . "/state");
	}
}

exit;