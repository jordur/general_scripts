#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with activities (BURROchratic approach)
# @Author: Arbol
# @Contributors:
########################################################

package activity;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use DateTime;
use Data::Dumper;

# -----------------------
# -- Global variables  --
# -----------------------

sub new
# Creation of new activity object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my ($sec,$min,$hour,$day,$month,$year,$week_day,$day_year,$daylight_savings) = localtime(time); 
	my %defaults = (
		start_time => DateTime->now(),
		stop_time => undef,
		name => "burocracia",
		dedication => 1, # Percentage of dedication
	);
	
	my $self = {%defaults,@_};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub set_parameter ()
# Function for setting parameters of activity objects
{
	my ($self, $param, $value) = @_;
	$self->{$param} = $value;
}

sub set_stop_time ()
# Function for setting parameters of activity objects
{
	my ($self) = @_;
	$self->{stop_time} = DateTime->now();
}

sub get_parameter ()
# Function for getting parameters of variants_info objects
{
	my ($self, $param) = @_;
	if (defined($self->{$param})){
		return $self->{$param};
	} else {
		return;
	}
}

sub get_ids ()
{
	my ($self) = @_;
	return keys (%{$self});
}

sub get_datetime($$){
	my ($date,$time) = @_;
	if (not defined($date) or not defined($time)){
		return undef;
	} else {
		my @date_array = split("-",$date);
		my @time_array = split(":",$time);
		my $epoch = DateTime->new(
			year	=> $date_array[2],
			month	=> $date_array[1],
			day		=> $date_array[0],
			hour	=> $time_array[0],
			minute	=> $time_array[1],
			second	=> $time_array[2],
		);
		return $epoch;
	}
}

sub get_activity_from_line
# Function that retrieves activity information from a tsv line
{
	my ($line) = @_;
	chomp($line);
	my @fields = split ("\t",$line);
	my $start = get_datetime($fields[2],$fields[3]);
	my $stop = get_datetime($fields[4],$fields[5]);
	my $self = activity->new(name => $fields[0], dedication => $fields[1],start_time => $start, stop_time => $stop);
	return $self;
}

sub get_activities_from_tsv
# Function that retrieves activities information from tsv files
{
	my ($class, $file_name) = @_;
	# Open file
	open (TSV, "<", $file_name) or die("ERROR: Couldn't open file $file_name!! $!\n");
	
	my %activities;
	my $num = 0;
	
	while (<TSV>){
		my $line = $_;
		if (substr($line,0,1) ne "#"){
			$activities{$num} = &get_activity_from_line($line);
			$num++;
		}
	}
	return \%activities;
}

sub transform_to_human{
	my ($epoch) = @_;
	if ($epoch ne ""){
		my $dt = DateTime->from_epoch( epoch => $epoch );
		return $dt->dmy . "\t". $dt->hms;
	} else {
		return "\t";
	}
}

sub put_to_line
# Function that exports an activity information to a tsv line
{
	my ($self) = @_;
	if (defined $self->{stop_time}){
		return $self->{name} . "\t" . $self->{dedication} . "\t" . &transform_to_human($self->{start_time}->epoch()) . "\t" . &transform_to_human($self->{stop_time}->epoch()) . "\n";
	} else {
		return $self->{name} . "\t" . $self->{dedication} . "\t" . &transform_to_human($self->{start_time}->epoch()) . "\t\t\n";
	}
}

sub print_to_file
# Function that prints an activity information to a tsv file
{
	my ($self, $file) = @_;
	my $line = $self->put_to_line();
	print $file $line;
}

sub put_activities_into_file
# Function that exports activities to a tsv file
{
	my ($class, $activities, $file_name) = @_;
	
	# Open file
	open (TSV, ">", $file_name) or die("ERROR: Couldn't open file $file_name!! $!\n");
	foreach (sort {$a <=> $b} keys (%{$activities})){
		my $activity = $$activities{$_};
		my $line = $activity->put_to_line();
		print TSV $line;
	}
}

1;