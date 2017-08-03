#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa
# @Contributors: Arbol
########################################################

package mapping_bwa;

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

sub mapping {
	my ($class, $essay, $jobs) = @_;
	print "Launching mapping_bwa jobs for essay $essay->{name}... ".localtime()."\n";
	my $list2;
	foreach my $sample_name (keys(%{$essay->get_parameter("samples")})){
		foreach my $replicate_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates"))})){
			my $list;
			foreach my $lane_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates",$replicate_name,"lanes"))})){
				my $lane_id = &aln_lane($essay,$jobs,$sample_name,$replicate_name,$lane_name);
				my $bam_id = sam_handle->manage_sam($essay,$lane_id,$sample_name,$replicate_name,$lane_name);
				$list = join_jobs($list,$bam_id);
			}
			my @files = @{$essay->get_bams(sample=>$sample_name,replicate=>$replicate_name,level=>"lane",type=>"PE")};
			my $merge_id = &merge_lanes (essay=>$essay,jobs=>$list,files=>\@files,sample=>$sample_name,replicate=>$replicate_name);
			my $hq_bam_id = &get_hq_bam($essay,$merge_id,$sample_name,$replicate_name);
			my $realign_id = &realign ($essay,$hq_bam_id,$sample_name,$replicate_name);
			my $recal_id =  &recalibrate ($essay, $realign_id,$sample_name,$replicate_name);
			$list2 = join_jobs($list2,$recal_id);
		}
	}
	#my @files = @{$essay->get_bams()};
	#my $merge_id = &merge_lanes (essay=>$essay,jobs=>$list2,files=>\@files);
	return $list2;
}

sub aln_lane {
	my ($essay,$jobs,$sample,$replicate,$lane) = @_;
	my @aln_id;
	my $sampe_id = "";
	if ($essay->get_lane_type($sample,$replicate,$lane) eq "PE"){
		my $left_id = &launch_bwa_aln ($essay,$jobs,$sample,$replicate,$lane,"R1","left");
		my $right_id = &launch_bwa_aln ($essay,$jobs,$sample,$replicate,$lane,"R2","right");
		my $aln_id = join_jobs ($left_id,$right_id);
		$sampe_id = &launch_bwa_sampe ($essay,$aln_id,$sample,$replicate,$lane);
	} else {
		die ("ERROR: Data type not yet implemmented!!");
	}
	return $sampe_id;
}

sub launch_bwa_sampe {
	my ($essay,$previous_jobs,$sample,$replicate,$lane) = @_;
	
	# Define involved information
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"bwa_sampe");
	my $reference = $essay->{references}{bwa_fasta};
	my $lane_R1 = $essay->get_reads($sample,$replicate,$lane,"R1");
	my $lane_R2 = $essay->get_reads($sample,$replicate,$lane,"R2");
	my $sai_left = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".left");
	my $sai_right = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".right");
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);
	my $sampe = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".sam");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "bwa sampe -s -a 250 -o 5 -P -n 5 -N 5 $reference $sai_left $sai_right $lane_R1 $lane_R2 -f $sampe";
	
	# Check if it's a panel or exome to determine the number of slots for sampe job.
	if (not defined($essay->{target_reference}{target}) or $essay->{target_reference}{target} =~ /Exoma/i){
		my $job = launch_essay_job ($job_file, '-pe smp 4', $previous_jobs, $string, $log_dir, $essay, $module);
	 } else {
	 	my $job = launch_essay_job ($job_file, '-pe smp 2', $previous_jobs, $string, $log_dir, $essay, $module);
	 }
	
	#my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub launch_bwa_aln {
	my ($essay,$previous_jobs,$sample,$replicate,$lane,$reads,$orientation) = @_;
	
	# Define involved files
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"alignment_".$orientation);
	my $reference = $essay->{references}{bwa_fasta};
	my $lane_file = $essay->get_reads($sample,$replicate,$lane,$reads);
	my $read_length = $essay->get_reads($sample,$replicate,$lane,"length_".$reads);
	my $mismatches = "-i ".int($read_length / 10);
	my $sai = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".".$orientation);
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);

	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "bwa aln $reference $lane_file -f $sai -o 20 -k 3 -l 30 -q 30 -t 4 $mismatches -e 50";
	my $job = launch_essay_job ($job_file, '-pe smp 4', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

1;