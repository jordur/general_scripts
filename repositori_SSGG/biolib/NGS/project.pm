#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with analysis folders
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

package project;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use NGS::essay;
use json_handle;
use Exporter;
use File::Which;
use File::Path qw(make_path remove_tree);

# -----------------------
# -- Global variables  --
# -----------------------

our @ISA = qw(Exporter);
our @EXPORT = qw ();

sub new
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;

	my %samples;
	my %references;
	my %essays;
		
	my %defaults = (
		name => "",
		path => "/",
		references => \%references,
		essays => \%essays,
		samples => \%samples,
	);
	
	my $self = {%defaults, @_};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub add_reference {
	my ($self, $name, $reference) = @_;
	$self->{references}{$name} = $reference;
}

sub add_essay {
	my ($self, $name, $essay) = @_;
	$self->{essays}{$name} = $essay;
}

sub add_sample {
	my ($self, $name, $sample) = @_;
	$self->{samples}{$name} = $sample;
}

sub get_project {
	my ($class, $json, $queue, $priority, $create_only_mode) = @_;
	my %data = &get_from_json($json);
	
	my %samples = %{$data{samples}};
	my %references = %{$data{references}};
	
	# Create essay hash
	my %essays;
	foreach my $essay_id (keys(%{$data{essays}})){
		my %target_reference;
		if (defined($data{essays}{$essay_id}{target_reference})){
			%target_reference = %{$data{essays}{$essay_id}{target_reference}};
		}
		my $pipeline = $data{essays}{$essay_id}{pipeline};
		
		# Modules to launch are defined in pipeline, but specific parameters may be overwritten
		# by those in project json
		my %modules = %{$data{essays}{$essay_id}{modules}};
		my $pipeline_json = which($data{essays}{$essay_id}{pipeline}.".json");
		if ($pipeline_json eq ""){
			die "ERROR: pipeline ".$data{essays}{$essay_id}{pipeline}." not available!!\n";
		}
		my %pipe = get_from_json($pipeline_json);
		# Initialize all modules in the pipeline
		foreach my $module (keys(%{$pipe{modules}})){
			if (not defined($modules{$module})){
				$modules{$module} = {};
			}
		}
		
		# Initialize modules state
		my $path = $data{path} . "/" . $essay_id;
		my $essay_tree;
		# Get essay state information, in case of existing
		if ( -e ($path . "/state.json")){
			my %existent_essay = get_from_json($path."/state.json");
			for my $current_module (keys(%modules)){
				# Update jobs states
				if (defined($existent_essay{modules}{$current_module}{state})){
					$modules{$current_module}{state} = $existent_essay{modules}{$current_module}{state};
				} else {
					$modules{$current_module}{state} = "pending";
				}
				for my $job (keys(%{$existent_essay{modules}{$current_module}{jobs}})){
					$modules{$current_module}{jobs}{$job} = $existent_essay{modules}{$current_module}{jobs}{$job};
				}
			}
			# Update tree information
			$essay_tree = $existent_essay{tree};
		} else {
			for (keys(%modules)){
				$modules{$_}{state} = "pending";
			}
			my %tree;
			$essay_tree = \%tree;
		}
		
		my @essay_samples_names = split(",",$data{essays}{$essay_id}{samples});
		my %essay_samples;
		foreach (@essay_samples_names){
			$essay_samples{$_} = $samples{$_};
		}
		my $essay = essay->new(
			name => $essay_id,
			path => $path,
			queue => $queue,
			priority => $priority,
			pipeline => $pipeline,
			seq_platform => "",
			references => \%references,
			target_reference => \%target_reference,
			modules => \%modules,
			samples => \%essay_samples,
			samples_order => \@essay_samples_names,
			tree => $essay_tree,
			create_only_mode => $create_only_mode,
			project_name => $data{name},
		);
		# Check species, pipeline and platform, and set corresponding values in essay object
		# Please note that following step redefines the hashtag "references" in the state.json
		$essay->checks(\%references);
		
		# Output log information related to repositories and pipeline version
		my $pipe_version = $essay->add_file("pipeline.version");
		my $pipe = $essay->get_parameter("pipeline");
		$essay->create_tree();
		print "INFO: Launching pipe_version.py -p $pipe > $pipe_version\n";
		my $output = `pipe_version.py -p $pipe > $pipe_version`;
				
		# Add essay to projects essays
		$essays{$essay_id} = $essay;
	}
	
	my $self = project->new(
		name => $data{name},
		path => $data{path},
		essays => \%essays,
		samples => \%samples,
	);
	return $self;
}

sub get_keys {
	my ($self) = shift;
	return (sort {$a cmp $b} keys %{$$self});
}

1;