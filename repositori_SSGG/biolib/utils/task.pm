#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with tasks (BURROchratic approach)
# @Author: Arbol
# @Contributors:
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;

# -----------------------
# -- Global variables  --
# -----------------------

package task;
sub new
# Creation of new tasks information object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my %defaults = (
		short => "",
		name => "",
		department => "",
		CC => "", # "cost center"
	);
	
	my $self = {%defaults,@_};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub set_parameter ()
# Function for setting parameters from samples of variants_info objects
{
	my ($self, $param, $value) = @_;
	$self->{$param} = $value;
}

sub get_parameter ()
# Function for getting parameters from samples of variants_info objects
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

sub get_task_from_line
# Function that retrieves task information from a tsv line
{
	my ($line) = @_;
	chomp($line);
	my @fields = split ("\t",$line);
	my $self = task->new(short => $fields[0],name => $fields[1], department => $fields[2], CC => $fields[3]);
	return $self;
}

sub get_tasks_from_tsv
# Function that retrieves tasks information from tsv files
{
	my ($class, $file_name) = @_;
	# Open file
	open (TSV, "<", $file_name) or die("ERROR: Couldn't open file $file_name!! $!\n");
	
	my %tasks;
	my $num = 0;
	
	while (<TSV>){
		my $line = $_;
		if (substr($line,0,1) ne "#"){
			$tasks{$num} = &get_task_from_line($line);
			$num++;
		}
	}
	return \%tasks;
}

sub put_to_line
# Function that exports a task information to a tsv line
{
	my ($self) = @_;
	return $self->{short} . "\t" . $self->{name} . "\t" . $self->{department} . "\t" . $self->{CC} . "\n";
}

sub print_to_file
# Function that prints an task information to a tsv file
{
	my ($self, $file) = @_;
	my $line = $self->put_to_line();
	print $file $line;
}

sub put_tasks_into_file
# Function that exports tasks to a tsv file
{
	my ($class, $tasks, $file_name) = @_;
	
	# Open file
	open (TSV, ">", $file_name) or die("ERROR: Couldn't open file $file_name!! $!\n");
	
	foreach (sort {$a <=> $b} keys (%{$tasks})){
		my $task = $$tasks{$_};
		my $line = $task->put_to_line();
		print TSV $line;
	}
}

1;