#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: Arbol
# @Contributors: 
########################################################

package primary_stats;

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
my $module = "primary_stats";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub stats
# Function for retrieving primary stats
{
	my ($class,$essay,$previous_jobs) = @_;
		
	print "Launching primary_stats jobs for essay $essay->{name}...".localtime()."\n";
	foreach my $sample (@{$essay->{samples_order}}){
		foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
			foreach my $lane (sort(keys(%{$essay->{samples}{$sample}{replicates}{$replicate}{lanes}}))){
				my $readsparam;
				my ($f3qual,$f5qual,$r1fastq,$r2fastq);
				if ($essay->{seq_platform} eq "SOLiD"){
					my $f3fasta = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F3_csfasta};
					$f3qual = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F3_qual};
					my $f5fasta = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F5_csfasta};
					$f5qual = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F5_qual};
					$readsparam = "-f3f $f3fasta -f3q $f3qual -f5f $f5fasta -f5q $f5qual";
					
				} elsif ($essay->{seq_platform} eq "Illumina"){
					$r1fastq = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{R1};
					$r2fastq = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{R2};
					$readsparam = "-r1 $r1fastq -r2 $r2fastq";
				}
				my $output = $essay->add_file("analysis","stats",$sample,$replicate,$lane,"rawdata_metrics.txt");
				
				# Define involved information
				my $job_file = $essay->add_file("jobs","stats",$sample,$replicate,$lane,"lane_checker");
				my $log_dir = $essay->mkdir("logs","stats",$sample,$replicate,$lane);
							
				# Create involved folders
				$essay->create_tree();
		
				# Launch jobs
				my $string = "lane_checker.pl $readsparam -o $output";
				my $job = launch_essay_job($job_file,"",$previous_jobs,$string,$log_dir,$essay,$module);
				
				# Define involved information
				$job_file = $essay->add_file("jobs","stats",$sample,$replicate,$lane,"primary_stats");
				
				if ($essay->{seq_platform} eq "SOLiD"){
					my $F3_stats_prefix = $essay->add_file("analysis","stats",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_F3");
					my $F5_stats_prefix = $essay->add_file("analysis","stats",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_F5");
					$string = "sg_primary_quality.sh $f3qual 1 no $F3_stats_prefix
					sg_primary_quality.sh $f5qual 1 no $F5_stats_prefix";
				} elsif ($essay->{seq_platform} eq "Illumina"){
					my $R1_stats = $essay->mkdir("analysis","stats",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_R1");
					my $R2_stats = $essay->mkdir("analysis","stats",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_R2");
					$string = "fastqc -o $R1_stats $r1fastq
					fastqc -o $R2_stats $r2fastq";
				}
							
				# Create involved folders
				$essay->create_tree();
		
				# Launch jobs
				$job = launch_essay_job($job_file,'',$previous_jobs,$string,$log_dir,$essay,$module);
					
			}
		}
	}
}

1;