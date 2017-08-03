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
	my ($job,$threads,$hold,$string,$log_dir,$essay,$module) = @_;
	my $col = '$3';
	my $task;
	
	# Update essay state file
	create_json($essay,$essay->{path} . "/state");
	
	# Check if JOB was already finished, and if so, job won't be relaunched
	my $job_state = $essay->GetJobState($module,$job);
	if (substr($job_state,0,8) eq "finished"){
		print "INFO: Task $job already executed and finished\n";
	} elsif ((substr($job_state,0,7) eq "running") or (substr($job_state,0,7) eq "procces")){
		
	} else {
		print "INFO: Launching task $job\n";
		chdir $log_dir or die "Can't move to $log_dir\n";
		my @names = split (/\//, $job);
		&create_job_file ($job, $string);
		$task =  Proc::Background->new($string);
		$essay->SetJobState($module,$job,"proccess=" . $task->{_pid});
	}
	return $task->{_pid};
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
	if (substr($essay->GetModuleState($module),0,8) ne "finished"){
		my $col = '$3';
		my $hold = "";
		
		# If current job is the first of the module, it has to wait for previous module jobs to be finished
		if (ref($jobs) ne "ARRAY" and $jobs eq ""){
			$hold = $essay->GetModuleHoldJobIDs($module);
			if ($hold ne ""){
				$hold = "-hold_jid " . $hold;
			}
		} elsif ( $#{$jobs} > -1) {
			foreach my $prev_job (@{$jobs}){
				my $current_hold = $essay->GetJobIDs($module,$prev_job);
				if ($current_hold ne ""){
					$hold = $hold . $current_hold . ",";
				}
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
			print "INFO: Job $job already executed and finished\n";
		} elsif ((substr($job_state,0,7) eq "running") or (substr($job_state,0,7) eq "procces")) {
			$job_id = substr($job_state,8,length($job_state));
			print "INFO: Not launching again job $job, since it's already running with ID $job_id\n";
		} else {
			print "INFO: Launching job $job (waiting for $hold jobs)\n";
			chdir $log_dir or die "Can't move to $log_dir\n";
			my @names = split (/\//, $job);
			&create_job_file ($job, $string);
			$job_id = `submit_job -q $essay->{queue} -p $essay->{priority} -N $names[$#names] $threads $hold $job.job | awk '{print $col}'`;
			chomp $job_id;
			$essay->SetJobState($module,$job,"running=" . $job_id);
		}
	}
	
	# Update essay state file
	create_json($essay,$essay->{path} . "/state");
	
	return $job_id;
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