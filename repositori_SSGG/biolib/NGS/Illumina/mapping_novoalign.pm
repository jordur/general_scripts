#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa
# @Contributors: Arbol
########################################################

package mapping_novoalign;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use NGS::essay;
use file_handle;
use sam_handle;
use bam_handle;
use submit;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------
my $module = "mapping";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub mapping{
	my ($class, $essay, $jobs_id) = @_;
	print "Launching mapping_novoalign jobs for essay $essay->{name}... ".localtime()."\n";
	my $list;
	foreach my $sample_name (keys(%{$essay->get_parameter("samples")})){
		foreach my $replicate_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates"))})){
			foreach my $lane_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates",$replicate_name,"lanes"))})){
				my $lane_id = &aln_lane($essay,$jobs_id,$sample_name,$replicate_name,$lane_name);
				my $bam_id = sam_handle->manage_sam($essay,$lane_id,$sample_name,$replicate_name,$lane_name);
				$list = join_jobs($list,$bam_id);
			}
		}
	}
	my @files = split (/\,/, $essay->get_lanes_bams());
	my $merge_id = &merge_lanes ($essay, $list, @files);
	my $realign_id = &realign ($essay, $merge_id);
	my $recal_id =  &recalibrate ($essay, $realign_id);
	my $hq_bam_id = &get_hq_bam($essay, $realign_id);
	return $hq_bam_id;
}

sub aln_lane {
	my ($essay,$jobs_id,$sample,$replicate,$lane) = @_;
	my $decomp_id = &decompress ($essay,$jobs_id,$sample,$replicate,$lane);
	return &launch_novoalign ($essay,$decomp_id,$sample,$replicate,$lane);
}

sub decompress {
	my ($essay, $previous_jobs,$sample,$replicate,$lane) = @_;
		
	# Define involved information
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"decompress");
	my $R1 = $essay->get_reads($sample,$replicate,$lane,"R1");
	my $R2 = $essay->get_reads($sample,$replicate,$lane,"R2");
	my $R1_decompres = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"{$sample}_{$replicate}_{$lane}_R1.fastq");
	my $R2_decompres = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"{$sample}_{$replicate}_{$lane}_R2.fastq");
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "gunzip -c $R1 > $R1_decompres
	gunzip -c $R2 > $R2_decompres";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub launch_novoalign {
	my ($essay,$previous_jobs,$sample,$replicate,$lane) = @_;
		
	# Define involved information
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"novoalign");
	my $novo_index = $essay->get_parameter("references","novoalign_index");
	my $R1_decompres = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"{$sample}_{$replicate}_{$lane}_R1.fastq");
	my $R2_decompres = $essay->add_file("trash","mapping",$sample,$replicate,$lane,"{$sample}_{$replicate}_{$lane}_R2.fastq");
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);
	my $sampe = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".sam");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "novoalign -d $novo_index -f $R1_decompres $R2_decompres -o SAM > $sampe";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
} 

1;