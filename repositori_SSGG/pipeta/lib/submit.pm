#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: Arbol, JM Rosa
# @Contributors: 
########################################################

package submit;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use file_handle;
use Exporter;
use Proc::Background;
use json_handle;

# -----------------------
# -- Global variables  --
# -----------------------

our @ISA = qw(Exporter);
our @EXPORT = qw (launch_job launch_task launch_essay_job join_jobs add_job_ids);

sub create_job_file {
	my ($file, $string) = @_;	
	my $output = &get_out_file_handle($file, 'job');
	print $output $string;
	return 1;
}

sub launch_task{
	my ($job,$threads,$jobs,$string,$log_dir,$essay,$module) = @_;
	
	# Variables definition
	my $col = '$3';
	my $task;
	my @hold_jobs;
	
	# Get jobid number (this ID is not an SGE id, that may change, it's an unmutable pipeta job id)
	my $id = $essay->SetJobID($module,$job);
	
	# If current job is the first of the module, it has to wait for previous module jobs to be finished
	if (ref($jobs) ne "ARRAY" and $jobs eq ""){
		@hold_jobs = $essay->GetModuleHoldJobIDs($module);
	} elsif ( $#{$jobs} > -1) {
		foreach my $prev_job (@{$jobs}){
			push(@hold_jobs,$essay->GetJobID($module,$prev_job));
		}
	} else {
		die "ERROR: Problem with jobs IDs that current job must wait for\n";
	}
	
	# Update essay state file
	create_json($essay,$essay->{path} . "/state");
	my $task_id = "";
	
	# Check if JOB was already finished, and if so, job won't be relaunched
	my $job_state = $essay->GetJobState($module,$job);
	if (substr($job_state,0,8) eq "finished"){
		$task_id = $essay->GetJobSGEID($module,$job);
		print "INFO: Task $job already executed and finished\n";
	} elsif (substr($job_state,0,7) eq "pending") {
		print "INFO: Launching task $job\n";
		chdir $log_dir or die "Can't move to $log_dir\n";
		my @names = split (/\//, $job);
		&create_job_file ($job, $string);
		if (not $essay->{'create_only_mode'}){
			$task =  Proc::Background->new($string);
			$essay->SetJobState($module,$job,"process");
			$task_id = $task->{_pid};
		} else {
			$essay->SetJobState($module,$job,"pending");
		}
	}
	$essay->SetJobID($module,$job,$task_id);
	$essay->SetJobSGEID($module,$job,$task_id);
	$essay->SetJobLog($module,$job,$log_dir);
	$essay->SetHoldJobs($module,$job,\@hold_jobs);
	return $id;
}

sub launch_job {
	my ($job, $threads, $hold, $queue, $priority, $string, $log_dir) = @_;
	my $col = '$3';
	chdir $log_dir or die "Can't move to $log_dir\n";
	my @names = split (/\//, $job);
	&create_job_file ($job, $string);
	my $job_id = `submit_job -q $queue -p $priority -N $names[$#names] $threads $hold $job.job | awk '{print $col}'`;
	chomp $job_id;	
	return $job_id;
}

sub launch_essay_job {
	my ($job,$threads,$jobs,$string,$log_dir,$essay,$module) = @_;
	my $job_id = "";
	
	# Check if MODULE was already finished, and if so, NO job from this module will be launched,
	# no matter what specific information for the individual job is available

	# Get jobid number (this ID is not an SGE id, that may change, it's an unmutable pipeta job id)
	my $id = $essay->SetJobID($module,$job);
	
	# Variables definitions
	my $col = '$3';
	my $hold = "";
	my @hold_jobs;
	
	# If current job is the first of the module, it has to wait for previous module jobs to be finished
	if (ref($jobs) ne "ARRAY" and $jobs eq ""){
		@hold_jobs = $essay->GetModuleHoldJobIDs($module);
		$hold = $essay->GetModuleHoldString($module);
	} elsif ( $#{$jobs} > -1) {
		foreach my $prev_job (@{$jobs}){
			my $current_hold = $essay->GetJobSGEID($module,$prev_job);
			if ($current_hold ne ""){
				$hold = $hold . $current_hold . ",";
			}
			push(@hold_jobs,$essay->GetJobID($module,$prev_job));
		}
		$hold =~ s/,$//;
		if ($hold ne ""){
			$hold = "-hold_jid " . $hold;
		}
	} else {
		die "ERROR: Problem with jobs IDs that current job must wait for\n";
	}
	
	# Check if JOB was already finished, and if so, job won't be relaunched
	my $job_state = $essay->GetJobState($module,$job);
	if (substr($job_state,0,8) eq "finished"){
		$job_id = $essay->GetJobSGEID($module,$job);
		print "INFO: Job $job already executed and finished\n";
	} elsif (substr($job_state,0,7) eq "pending") {
		print "INFO: Launching job $job (waiting for $hold jobs)\n";
		chdir $log_dir or die "Can't move to $log_dir\n";
		my @names = split (/\//, $job);
		&create_job_file ($job, $string);
		if (not $essay->{'create_only_mode'}){
			# TO_REMOVE: Temporary changes for allowing multi-thread jobs (and also jobs that shouldn't be suspended by threshold) to be addressed to higher priority queues
			if ($names[$#names]=~/ensembl/ or $names[$#names]=~/merge/){
				$job_id = `submit_job -q vaishia.q -p $essay->{priority} -N $names[$#names] $threads $hold $job.job | awk '{print $col}'`;
			} else {
				$job_id = `submit_job -q $essay->{queue} -p $essay->{priority} -N $names[$#names] $threads $hold $job.job | awk '{print $col}'`;
			}
			chomp $job_id;
			$essay->SetJobState($module,$job,"running");
		} else {
			$essay->SetJobState($module,$job,"pending");
		}
	} else {
		$job_id = $essay->GetJobSGEID($module,$job);
		print "INFO: Not launching again job $job, since it's already active with ID $job_id\n";
	}
	$essay->SetJobSGEID($module,$job,$job_id);
	$essay->SetJobLog($module,$job,$log_dir);
	$essay->SetHoldJobs($module,$job,\@hold_jobs);
	
	# Update essay state file
	create_json($essay,$essay->{path} . "/state");
	return $id;
}

sub join_jobs{
	my ($jobs_array,$jobs_to_add) = @_;
	foreach (@{$jobs_to_add}){
		push(@{$jobs_array},$_);
	}
	return $jobs_array;
}

sub add_job_ids{
	my ($jobs_to_wait,$job_list) = @_;
	foreach (@{$job_list}){
		push(@{$jobs_to_wait},$_);
	}
	return $jobs_to_wait;
}

1;