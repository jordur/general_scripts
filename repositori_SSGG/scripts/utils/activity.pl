#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use POSIX;
use Data::Dumper;
use lib '/share/apps/scripts/global/SSGG_libraries';
use utils::activity;
use utils::task;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $activities_names;
my $start = 0;
my $finish = 0;
my $date = "";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"a=s"	=> \$activities_names,
				"s"		=> \$start,
				"f"		=> \$finish,
				"d=s"		=> \$date,
);

if (not defined($activities_names) and $finish == 0) {
	die "\n=================================================================
$basename: Script for registering dedications to different activities v1.1\n		  
	  Options:
	-a activ1,activ2  activity or colon separated list of activities that are going to be carried out
	-s                start time counter. Use this action at the beginning of the day, in case that tasks which you will dedicate to, remain unknown
	-f                finish all open activities. Use this action at the end of the day, at lunch time, etc.
	-d date           date whose tasks are going to be modified from (default is current date). Format: yyyy-mm-dd
	  Output:
	<USER_HOME>/repercusiones/YYYYMM.txt  File containing all performed activities during a month, in following format:
	activity_shortening <TAB> dedication <TAB> start_date(yyyy-mm-dd) <TAB> start_time(hh:mm:ss) <TAB> stop_date <TAB> stop_time_in_seconds
	For a complete list of available tasks, check /share/gluster/burocracia/tareas.txt
\n".localtime()."\n=================================================================\n\n";
}

&main ();
print  "\nINFO: Activity added successfully... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	my $path = `pwd`;
	chomp $path;
	
	my $dedications_folder = $ENV{"HOME"} . "/repercusiones";
	unless(-d $dedications_folder){
		mkdir $dedications_folder or die ("ERROR: Impossible to create folder $dedications_folder!!"); 
	}
	my $yyyymm = strftime("%Y%m", localtime);
	my $month_file_name = $dedications_folder . "/" . $yyyymm .".txt";
	
	my $tasks_file_name = "/share/gluster/burocracia/tareas.txt";
	
	# Open file
	open (my $day_file, ">>", $month_file_name) or die("ERROR: Couldn't open file $month_file_name!! $!\n");
	
	# Load tasks information
	my $tasks = task->get_tasks_from_tsv($tasks_file_name);
	
	# Load all month activities
	my $activities;
	if ( -e $month_file_name ){
		$activities = activity->get_activities_from_tsv($month_file_name);
	}
	
	if ($start){
		# Start a new activity with "empty" name
		
	}
	if ($finish){
		# Close all open activities
		foreach (keys (%{$activities})){
			my $activity_id = $_;
			# Set stop times for unconcluded activities. Names will be later added to "open"
			# activities (those with an empty name)
			if ($$activities{$activity_id}->get_parameter("name") ne ""){
				if (not defined($$activities{$activity_id}->get_parameter("stop_time"))){
					$$activities{$activity_id}->set_stop_time();
				}
			} else {
				my $start_time = $$activities{$activity_id}->get_parameter("start_time");
				print "WARNING: Open activity started at $start_time without task name assigned. No stop times will be set for this activity!";
			}
		}
	} else {
		# Get all new activity names, number of new activities and dedication
		my @activs = split (',',$activities_names);
		my $activs_number = $#activs + 1;
		my $dedication = 1 / $activs_number;
		
		# Check if all activities in list are already in tasks list
		foreach (@activs){
			my $activity_name = $_;
			my $found = 0;
			foreach (keys(%{$tasks})){
				if ($$tasks{$_}->get_parameter("short") eq $activity_name){
					$found = 1;
				}
			}
			if (not $found){
				print "ERROR: Nombre de actividad desconocido: $activity_name\n\n";
				print "Nombres permitidos de tareas (ver $tasks_file_name):\n";
				foreach (keys %{$tasks}){
					my $short = $$tasks{$_}->get_parameter("short");
					my $name = $$tasks{$_}->get_parameter("name");
					print "$short : $name\n";
				}
				exit;
			}
		}
		
		# Store all new activities and close (set stop times) all open activities
		foreach (@activs){
			# Create activity object
			my $activity = activity->new(name => "$_",dedication => $dedication);
			my $activity_name = $activity->get_parameter("name");
			
			# An activity stop time is set when a new activity is started, 
			# or when a finish signal is sended,
			# or when a defined activity is manually included just after the activity
			
			# Get unconcluded activities and last activity_id
			my $last_activity_id = -1;
			foreach (keys (%{$activities})){
				my $activity_id = $_;
				# Set stop times for unconcluded activities
				if (not defined($$activities{$activity_id}->get_parameter("stop_time"))){
					$$activities{$activity_id}->set_stop_time();
				}
				if ($last_activity_id < $activity_id){
					$last_activity_id = $activity_id;
				}
			}
			# Add activity to all day activities
			$$activities{$last_activity_id + 1} = $activity;
		}
	}
	activity->put_activities_into_file($activities,$month_file_name);
}

exit;