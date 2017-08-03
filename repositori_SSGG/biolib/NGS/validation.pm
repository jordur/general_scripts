#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: Arbol
# @Contributors: 
########################################################

package validation;

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
my $module = "validation";
our @ISA = qw(Exporter);
our @EXPORT = qw ();

sub validate {
	my ($class, $essay, $id) = @_;
	print "Launching validation jobs for essay $essay->{name}... ".localtime()."\n";
	my $job_list;
	foreach my $validation (keys(%{$essay->{modules}{$module}{validations}})){
		# Define involved information
		my $log_dir = $essay->mkdir("logs","validation",$validation);
		my $variants = $essay->add_file("analysis","variants",$essay->{name}."_collect_not_annotated.vcf");
		my $panel_name = $essay->{target_reference}{name};
		my $samples = $essay->{modules}{$module}{validations}{$validation}{samples};
		my $hapmap_cell_line = $essay->{modules}{$module}{validations}{$validation}{HapMap_cell_line};
		my $output_path = $essay->mkdir("analysis","validation",$validation);
		my $tmp_folder = $essay->mkdir("trash","validation",$validation);
		# Launch jobs
		my $job_file1 = $essay->add_file("jobs","validation",$validation,"concordance");
		$essay->create_tree();
		my $string = "sg_concordance_studies.sh -v $variants -s $samples -d $panel_name -h $hapmap_cell_line -t $tmp_folder -o $output_path";
		my $job = launch_essay_job ($job_file1, "", $id, $string, $log_dir, $essay, $module);
		$job_list = join_jobs($job_list,[$job_file1]);
	}
	return $job_list;
}

1;