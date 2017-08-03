#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bowtie2
# @Author: Guillermo Marco
# @Contributors: Arbol
########################################################

package mapping_bowtie2;

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
	print "Launching mapping_bowtie2 jobs for essay $essay->{name}... ".localtime()."\n";
	my $list2;
	foreach my $sample_name (keys(%{$essay->get_parameter("samples")})){
		foreach my $replicate_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates"))})){
			my $list;
			foreach my $lane_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates",$replicate_name,"lanes"))})){
				my $lane_id = &bowtie2_aln($essay,$jobs,$sample_name,$replicate_name,$lane_name);
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

sub bowtie2_aln {
	my ($essay,$previous_jobs,$sample,$replicate,$lane) = @_;
	my ($job_file, $string, $lane_R1, $lane_R2);
	my $analysis_name = $essay->{'target_reference'}{'name'};
	if ($essay->get_lane_type($sample,$replicate,$lane) eq "PE"){
		# Define involved files
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,$lane,"bowtie2_alignment");
		my $reference = $essay->{references}{bowtie2_reference};
		
		#Check trimming if module is defined and if there's then use trimmed fastq.
		if ($essay->exists_module("trimming")){
			$lane_R1 = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R1_5_3_prime_paired_trimmed");
			$lane_R2 = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R2_5_3_prime_paired_trimmed")
		}
		else {
		$lane_R1 = $essay->get_reads($sample,$replicate,$lane,"R1");
		$lane_R2 = $essay->get_reads($sample,$replicate,$lane,"R2");			
		}

		my $trash_mapping = $essay->add_file("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".sam");
		my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate,$lane);

	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "bowtie2 -p 4 --very-sensitive -x $reference -1 $lane_R1 -2 $lane_R2 -S $trash_mapping";
	
	#BRCAs nextera
	if ($analysis_name eq 'brcas_nextera'){
		$string = "bowtie2 -p 4 --very-sensitive -I 250 -X 1000 -x $reference -1 $lane_R1 -2 $lane_R2 -S $trash_mapping";
	}
	
	my $job = launch_essay_job ($job_file, "-pe smp 4", $previous_jobs, $string, $log_dir, $essay, $module);
		
		
	} else {
		die ("ERROR: Data type not yet implemmented!!");
	}
	return [$job_file];
}

1;