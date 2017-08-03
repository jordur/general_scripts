#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map data
# @Author: JM Rosa, Arbol
# @Contributors: 
########################################################

# package that handles bam files once they have been collected and merged
package bam_handle;

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use File::Which;
use NGS::essay;
use file_handle;
use submit;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------

my $module = "mapping";
our @ISA = qw(Exporter);
our @EXPORT = qw(merge_lanes realign recalibrate get_hq_bam);

sub get_hq_bam {
	my ($essay,$id,$sample,$replicate,$type) = @_;
	my ($remove_q1_id,$pcr_id,$fixmate_id);
	
	my $current_capture_system = $essay->{samples}{$sample}{replicates}{$replicate}{capture_system};
	# NO FixMate will be carried out until problems with PicardTools are solved
	# It has been observed that FixMate retrieves UNPAIRED reads (with different IDs) and set
	# them together!!!!
	if (defined($type)){
		$remove_q1_id = &remove_Q1 ($essay,$id,$sample,$replicate,$type);
		$pcr_id = &remove_PCR ($essay,$remove_q1_id,$sample,$replicate,$type);
		#$fixmate_id = &fix_mate_info ($essay,$pcr_id,$sample,$replicate,$type);
	} else {
		
		#If capture system is Multiplicom we don't remove PCR duplicates.
		if ($current_capture_system =~ /Multiplicom/i || $current_capture_system =~ /Haloplex/i){
			$remove_q1_id = &remove_Q1 ($essay,$id,$sample,$replicate);
			return $remove_q1_id;
		}
		#If capture system is different than Multiplicom
		else {
			$remove_q1_id = &remove_Q1 ($essay,$id,$sample,$replicate);
			$pcr_id = &remove_PCR ($essay,$remove_q1_id,$sample,$replicate);	
		}
		#$fixmate_id = &fix_mate_info ($essay,$pcr_id,$sample,$replicate);
	}
	return $pcr_id;
}

sub fix_mate_info {
	my ($essay,$previous_jobs,$sample,$replicate,$type) = @_;
		
	# Define involved information
	my ($job_file,$input, $output);
	if ($essay->is_solid() eq 'true'){
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"fix_mate_info_".$type);
		if ($type eq "F3"){
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR_F3.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR_fixmate_F3.bam");
		} else {
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR_fixmate.bam");
		}
	}else {
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"fix_mate_info");
		$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR.bam");
		$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_PCR_fixmate.bam");
	}

	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	# Info: Parameter MAX_RECORDS_IN_RAM is set to 400000 (former value was 100000), because otherwise this step raises an error (too many open files) for big datasets (Exomes)
	my $string = "FixMateInformation.jar I=$input O=$output TMP_DIR=/scratch SO=coordinate CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=400000";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub remove_PCR {
	my ($essay,$previous_jobs,$sample,$replicate,$type) = @_;
		
	# Define involved information
	my $job_file;
	my ($input, $output);
	if ($essay->is_solid() eq 'true'){
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"remove_PCR_".$type);
		if ($type eq "F3"){
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_F3.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_F3.bam");
		} else {
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR.bam");
		}
	}else {
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"remove_PCR");
		$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");
		$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR.bam");
	}	 

	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	# Info: Parameter MAX_RECORDS_IN_RAM is set to 400000 (former value was 100000), because otherwise this step raises an error (too many open files) for big datasets (Exomes)
	my $string = "MarkDuplicates.jar I=$input O=$output M=stat TMP_DIR=/scratch REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE ASSUME_SORTED=FALSE MAX_RECORDS_IN_RAM=400000";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub remove_Q1 {
	my ($essay,$previous_jobs,$sample,$replicate,$type) = @_;
		
	# Define involved information
	my $job_file;
	my ($input, $output);
	if ($essay->is_solid() eq 'true'){
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"remove_Q1_".$type);
		if ($type eq "F3"){
			$input = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_F3.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_F3.bam");			
		} else {
			$input = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge.bam");
			$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");			
		}
	}else {
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"remove_Q1");
		$input = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge.bam");
		$output = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");
	}	 

	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "samtools view -bh -q 1 $input  > $output
	samtools index $output";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub recalibrate {
	my ($essay,$previous_jobs,$sample,$replicate,$type) = @_;
	
	my $realign="";
	my $current_capture_system = $essay->{samples}{$sample}{replicates}{$replicate}{capture_system};
		
	# Define involved information
	my $job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"recalibration");
	my $reference = $essay->{references}{fasta};
	my $vcf = $essay->{references}{vcf};
	my $target = "";
	if (defined($essay->{'target_reference'}{'target'})){
		$target = "-L ".$essay->{'target_reference'}{'target'};
	}
	
	#If capture system is Multiplicom.
	if ($current_capture_system =~ /Multiplicom/i || $current_capture_system =~ /Haloplex/i){
		$realign = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_realign.bam");
	}
	
	else{
		$realign = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign.bam");
	}
	
	#my $realign = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign.bam");
	my $recalibration = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_recaldata.grp");
	my $recal = $essay->add_file("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam");
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	# TODO: Remove the -rf BadCigar as soon as possible
	my $string = "gatk_toolkit -T BaseRecalibrator -rf BadCigar -I $realign -R $reference -knownSites $vcf --disable_indel_quals -o $recalibration $target
	gatk_toolkit -T PrintReads -R $reference -I $realign -o $recal -BQSR $recalibration";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub realign {
	my ($essay,$previous_jobs,$sample,$replicate,$type) = @_;
		
	# Define involved information
	my $job_file;
	my ($input,$realign,$intervals);
	my $current_capture_system = $essay->{samples}{$sample}{replicates}{$replicate}{capture_system};
	
	if (($essay->is_solid() eq 'true') and ($type eq "F3")){
		# For SOLiD - F3
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"realign_".$type);
		$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_F3.bam");
		$realign = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign_F3.bam");
		$intervals = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_F3.intervals");
	} elsif($essay->is_solid() eq 'true') {
		# For SOLiD - PE
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"realign");
		$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR.bam");
		$realign = $essay->add_file("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam");
		$intervals = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate.".intervals");
	} else {
		# For Illumina
		$job_file = $essay->add_file("jobs","mapping",$sample,$replicate,"realign");
		
		#Check if Multiplicom is current_capture_system: Remove Q1
		if ($current_capture_system =~ /Multiplicom/i){
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");
			$realign = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_realign.bam");
		}
		
		#If we're not using Multiplicom: Remove Q1+PCR
		else{
			$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR.bam");
			$realign = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign.bam");	
		}
		#$input = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR.bam");
		#$realign = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign.bam");
		$intervals = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate.".intervals");
	}
	my $reference = $essay->{'references'}{'fasta'};
	my $vcf = $essay->{'references'}{'vcf'};
	my $target = "";
	if (defined($essay->{'target_reference'}{'target'})){
		$target = "-L ".$essay->{'target_reference'}{'target'};
	}
	my $log_dir = $essay->mkdir("logs","mapping",$sample,$replicate);

	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "gatk_toolkit -T RealignerTargetCreator -I $input -R $reference -known:dbsnp,vcf $vcf -o $intervals $target
	gatk_toolkit -T IndelRealigner -targetIntervals $intervals -R $reference -maxReads 1000000 -filterMBQ -I $input -o $realign";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub merge_lanes {
	# This function merges bam lanes coming from types ($type) "replicate" or "essay"
	my %defaults = (essay=>undef,jobs=>undef,files=>undef,sample=>undef,replicate=>undef,type=>undef);
	my $parameters = { %defaults, @_};
	my $essay = $$parameters{essay};
	my $jobs = $$parameters{jobs};
	my $files = $$parameters{files};
	my $merge_sam_jar_path = which('MergeSamFiles.jar');
	
	# Define involved information
	my ($job_file,$merge,$log_dir);
	my $files_txt = join (' I=', @{$files});
	$files_txt = 'I=' . $files_txt;
	if (defined($$parameters{sample}) and defined($$parameters{replicate})) {
		if (($essay->is_solid() eq 'true') and $$parameters{type} eq "F3"){
			$job_file = $essay->add_file("jobs","mapping",$$parameters{sample},$$parameters{replicate},"merge_F3");
			$merge = $essay->add_file("trash","mapping",$$parameters{sample},$$parameters{replicate},$$parameters{sample}."_".$$parameters{replicate}."_merge_F3.bam");
		} else {
			$job_file = $essay->add_file("jobs","mapping",$$parameters{sample},$$parameters{replicate},"merge");
			$merge = $essay->add_file("trash","mapping",$$parameters{sample},$$parameters{replicate},$$parameters{sample}."_".$$parameters{replicate}."_merge.bam");
		}
		$log_dir = $essay->mkdir("logs","mapping",$$parameters{sample},$$parameters{replicate});
	} else {
		if (($essay->is_solid() eq 'true') and $$parameters{type} eq "F3"){
			$job_file = $essay->add_file("jobs","mapping","bam","merge_F3");
			$merge = $essay->add_file("analysis","mapping","bam",$essay->{name}."_merge_F3.bam");
		} else {
			$job_file = $essay->add_file("jobs","mapping","bam","merge");
			$merge = $essay->add_file("analysis","mapping","bam",$essay->{name}."_merge.bam");
		}
		$log_dir = $essay->mkdir("logs","mapping","bam");
	}
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	# Info: Parameter MAX_RECORDS_IN_RAM is set to 400000 (former value was 100000), because otherwise this step raises an error (too many open files) for big datasets (Exomes)
	my $string = "java -d64 -Xmx8g -jar $merge_sam_jar_path $files_txt O=$merge TMP_DIR=/scratch SO=coordinate CREATE_INDEX=TRUE ASSUME_SORTED=FALSE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=400000";
	my $job = launch_essay_job ($job_file, "-pe smp 2", $jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

1;