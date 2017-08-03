#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa
# @Contributors: Arbol
########################################################

package annotation;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use file_handle;
use submit;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------
my $module = "annotation";
our @ISA = qw(Exporter);
our @EXPORT = qw ();

sub annotate {
	my ($class, $essay, $id) = @_;
	print "Launching annotation jobs for essay $essay->{name}... ".localtime()."\n";
	my $anno_id = &ensembl($essay, $id);
	return $anno_id; #&collecting($essay, $anno_id);
}

sub collecting {
	my ($essay, $previous_jobs) = @_;
		
	# Define involved information
	my $job_file = $essay->add_file("jobs","annotation","final_collecting");
	my $annotation = $essay->add_file("analysis","annotation",$essay->{name}."_annotation");
	my $ensembl = $essay->add_file("trash","annotation",$essay->{name}."_ensembl");
	my $filtering = "";
	if ($essay->get_parameter("modules","annotation","parameters","filter") eq "yes"){
		$filtering = "--filtering";
	}
	my $log_dir = $essay->mkdir("logs","annotation");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "sg_parsing_freqs_pops.pl -i $ensembl -o $annotation $filtering";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub ensembl {
	my ($essay, $previous_jobs) = @_;
		
	# Define involved information
	my $job_file = $essay->add_file("jobs","annotation","ensembl");
	my $annotation = $essay->add_file("analysis","annotation",$essay->{name}."_annotation.vcf");
	my $filtered = $essay->add_file("analysis","variants",$essay->{name}."_collect.vcf");
	my $log_dir = $essay->mkdir("logs","annotation");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "vep_launcher.pl -t 6 -i $filtered -o $annotation";
	my $job = launch_essay_job ($job_file, "-pe smp 4", $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

1;