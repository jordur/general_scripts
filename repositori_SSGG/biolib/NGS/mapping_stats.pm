#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

package mapping_stats;

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
my $module = "mapping_stats";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub stats
# Function for retrieving coverage and mapping stats
{
	my ($class,$essay,$previous_jobs) = @_;
	
	print "Launching mapping_stats jobs for essay $essay->{name}...".localtime()."\n";
	my @wait_jobs;
	foreach my $sample (@{$essay->{samples_order}}){
		foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
			
			# Define involved information
			my $job_file = $essay->add_file("jobs","stats",$sample,$replicate,"stats_mapping");
			my $bam = $essay->add_file("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam");
			my $coverages = $essay->get_parameter("modules","mapping_stats","parameters","bamstats","depth_thresholds","graph_coverages_values");
			my $output = $essay->add_file("analysis","stats",$sample,$replicate,"mapping");
			my $log_dir = $essay->mkdir("logs","stats",$sample,$replicate);
			my $target = "";
			if (defined($essay->get_parameter("target_reference","target"))){
				$target = $essay->get_parameter("target_reference","target");
			}
			my $queue = $essay->{queue};
			my $priority = $essay->{priority};
						
			# Create involved folders
			$essay->create_tree();
	
			# Launch jobs
			my $hold_jobs = $essay->GetModuleHoldString($module);
			$hold_jobs =~ s/^-hold_jid//;
			my $string = "sg_mapping_stats.sh $bam $coverages 1 16 $output $queue $target $priority $hold_jobs";
			my $job = launch_task($job_file,'',$previous_jobs,$string,$log_dir,$essay,$module);
			
			# Mapping performance
			my @rawdata_metrics_reports;
			foreach my $lane (sort(keys(%{$essay->{samples}{$sample}{replicates}{$replicate}{lanes}}))){
				push (@rawdata_metrics_reports,$essay->get_path("analysis","stats",$sample,$replicate,$lane,"rawdata_metrics.txt"));
			}
			my $rawmetrics = join(",",@rawdata_metrics_reports);
			my $merge = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge.bam");
			my $pcr = $essay->get_path("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam");
			my $q1 = $essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1.bam");
			my $tmp_folder = $essay->mkdir("trash","stats",$sample,$replicate);
			$job_file = $essay->add_file("jobs","stats",$sample,$replicate,"mapping_performance");
			$log_dir = $essay->mkdir("logs","stats",$sample,$replicate);
			$output = $essay->add_file("analysis","stats",$sample,$replicate,"mapping_performance.txt");
			if ($target ne ""){
				$target = "-t ".$target;
			}
			$string = "mapping_checker.pl -r $rawmetrics -m $merge -q1 $q1 -pcr $pcr $target -o $output -f $tmp_folder";
			
			# Create involved folders
			$essay->create_tree();
	
			$job = launch_essay_job ($job_file,"",$previous_jobs,$string,$log_dir,$essay,$module);
			
			# Non covered bases stats
			my $nc_stats_folder = $essay->add_file("analysis","stats",$sample,$replicate,"NC_stats");
			my $panel_name = $essay->get_parameter("target_reference","name");
			$log_dir = $essay->mkdir("logs","stats",$sample,$replicate);
#			$job_file = $essay->add_file("jobs","stats",$sample,$replicate,"NC_stats");
			$essay->create_tree();
			
			my $nc_depths = $essay->get_parameter("modules","mapping_stats","parameters","NC_stats","depth_thresholds");
			my @depths = split (',', $nc_depths);

			#For each depth read from JSON config
			foreach my $depth (@depths){
				$job_file = $essay->add_file("jobs","stats",$sample,$replicate,"NC_stats_$depth");
				$string = "sg_panel_non_covered_stats.pl -i $pcr -prefix NC -d $depth -threshold 100 -o $nc_stats_folder -p $panel_name\n";
				$job = launch_essay_job ($job_file,"-pe smp 2",$previous_jobs,$string,$log_dir,$essay,$module);
			} 
#			$job = launch_essay_job ($job_file,"",$previous_jobs,$string,$log_dir,$essay,$module);
		}
	}
	return \@wait_jobs;
}

1;