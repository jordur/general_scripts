#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bowtie2
# @Author: Guillermo Marco
# @Contributors: Arbol
########################################################

package trimming;

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
my $module = "trimming";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub trimming {
	my ($class, $essay, $jobs) = @_;
	print "Launching trimming jobs for essay $essay->{name}... ".localtime()."\n";
	my $list;
	foreach my $sample_name (keys(%{$essay->get_parameter("samples")})){
		foreach my $replicate_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates"))})){
			
			foreach my $lane_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates",$replicate_name,"lanes"))})){
				my $lane_id = &trimmomatic($essay,$jobs,$sample_name,$replicate_name,$lane_name);
				$list = join_jobs($list,$lane_id);
			}
		}
	}
	return $list;
}

sub trimmomatic {
	my ($essay,$previous_jobs,$sample,$replicate,$lane) = @_;
	my ($job_file, $string);
	if ($essay->get_lane_type($sample,$replicate,$lane) eq "PE"){
		# Define involved files
		$job_file = $essay->add_file("jobs","trimming",$sample,$replicate,$lane,"trimming");
		
		#References
		my $five_prime_primers_reference = $essay->{target_reference}{five_prime_primers};
		my $three_prime_primers_reference = $essay->{target_reference}{three_prime_primers};
		
		#Rawdata input
		my $lane_R1 = $essay->get_reads($sample,$replicate,$lane,"R1");
		my $lane_R2 = $essay->get_reads($sample,$replicate,$lane,"R2");
		
		#5 prime trimming trash
		my $R1_trash_trimming_5_prime_paired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R1_5_prime_paired_trimmed");
		my $R2_trash_trimming_5_prime_paired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R2_5_prime_paired_trimmed");
		my $R1_trash_trimming_5_prime_unpaired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R1_5_prime_unpaired_trimmed");
		my $R2_trash_trimming_5_prime_unpaired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R2_5_prime_unpaired_trimmed");
		
		#3 prime trimming trash
		my $R1_trash_trimming_5_3_prime_paired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R1_5_3_prime_paired_trimmed");
		my $R2_trash_trimming_5_3_prime_paired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R2_5_3_prime_paired_trimmed");
		my $R1_trash_trimming_5_3_prime_unpaired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R1_5_3_prime_unpaired_trimmed");
		my $R2_trash_trimming_5_3_prime_unpaired = $essay->add_file("trash","trimming",$sample,$replicate,$lane,"R2_5_3_prime_unpaired_trimmed");
		
		$string = "java org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33 $lane_R1 $lane_R2 $R1_trash_trimming_5_prime_paired $R1_trash_trimming_5_prime_unpaired $R2_trash_trimming_5_prime_paired $R2_trash_trimming_5_prime_unpaired ILLUMINACLIP:$five_prime_primers_reference:2:5:5 MINLEN:30\n";		
		$string .= "java org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33 $R1_trash_trimming_5_prime_paired $R2_trash_trimming_5_prime_paired $R1_trash_trimming_5_3_prime_paired $R1_trash_trimming_5_3_prime_unpaired $R2_trash_trimming_5_3_prime_paired $R2_trash_trimming_5_3_prime_unpaired ILLUMINACLIP:$three_prime_primers_reference:2:5:5 MINLEN:30";
		
		my $log_dir = $essay->mkdir("logs","trimming",$sample,$replicate,$lane);

	#Create involved folders
	$essay->create_tree();
	
	#Launch jobs
	#my $string = "bowtie2 -p 4 --very-sensitive -x $reference -1 $lane_R1 -2 $lane_R2 -S $trash_mapping";
	my $job = launch_essay_job ($job_file, "-pe smp 4", $previous_jobs, $string, $log_dir, $essay, $module);
		
		
	} else {
		die ("ERROR: Data type not yet implemmented!!");
	}
	return [$job_file];
}

1;