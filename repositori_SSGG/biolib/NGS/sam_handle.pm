#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa, Arbol
# @Contributors: 
########################################################

# package that manages sam and bam files at individual lane level
package sam_handle;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use NGS::essay;
use file_handle;
use submit;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------
my $module = "mapping";
our @ISA = qw(Exporter);
our @EXPORT = qw(manage_sam);

sub manage_sam {
	my ($class,$essay,$id,$sample,$replicate,$lane,$type) = @_;
	my ($rg_id, @files, @info);
	my $clean_id;
	if ($essay->is_solid() eq 'true') {
		$clean_id = $id;
		$rg_id = &replace_read_group ($essay,$clean_id,$sample,$replicate,$lane,$type);
	} else {
		$clean_id = &remove_reads_sam ($essay,$id,$sample,$replicate,$lane);
		$rg_id = &replace_read_group ($essay,$clean_id,$sample,$replicate,$lane);
	}
	return $rg_id;
}

sub replace_read_group {
	my ($essay,$previous_jobs,$sample,$replicate,$lane,$type) = @_;
	
	# Define involved information
	# In case of SOLiD analysis, also F3 reads will be separatelly mapped
	# This parameter comes in the $type argument
	my ($input,$output,$rgsm,$rgpu,$rgid,$job_file);
	if ($essay->is_solid() eq 'true'){
		if ($type eq "PE"){
			$input = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"output","pairing","F3-F5-Paired.bam");
			$output =  $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".readgroup.bam");
		} else {
			# Look for the F3 bam file name
			opendir(DIR, $essay->mkdir("trash","mapping",$sample,$replicate,$lane) . "/output/F3/maToBam");
			my @files= grep {/\.csfasta\.ma\.bam$/} readdir DIR;
			if (defined($files[0]) and $files[0] ne ""){
				$input = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"output","F3","maToBam",$files[0]);
			} else {
				die "ERROR: BioScope bam files haven't been found in ".$essay->get_path("trash","mapping",$sample,$replicate,$lane) . "/output/F3/maToBam\n";
			}
			$output =  $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_".$type.".readgroup.bam");
		}
		$rgsm = $sample.".".$replicate.".".$type;
		$rgpu = $sample.".".$replicate.".".$lane.".".$type;
		$rgid = $sample.".".$replicate.".".$lane.".".$type;
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"replace_read_group_$type");
	} else {
		$input = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_filtered.sam");
		$output = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".readgroup.bam");
		$rgsm = $sample . "." . $replicate;
		$rgpu = $sample . "." . $replicate . "." . $lane;
		$rgid = $sample . "." . $replicate . "." . $lane;
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"replace_read_group");
	}
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);

	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "AddOrReplaceReadGroups.jar I=$input O=$output TMP_DIR=/scratch RGLB=SG RGPL=$essay->{seq_platform} RGSM=$rgsm RGPU=$rgpu RGID=$rgid CREATE_INDEX=TRUE";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub remove_reads_sam {
	my ($essay,$previous_jobs,$sample,$replicate,$lane) = @_;
	
	# Define involved information
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"remove_reads_sam");
	my $sampe = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".sam");
	my $filtered = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane);
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "sg_remove_reads_sam.pl -i $sampe  -o $filtered";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

1;