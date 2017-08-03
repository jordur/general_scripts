#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to define
# @Desc: NGS data analysis involved paths and files
# @Author: Arbol & JM Rosa
# @Contributors:
########################################################

package essay;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use file_handle;
use json_handle;
use Exporter;
use JSON -support_by_pp;
use File::Which;
use File::Path qw(make_path remove_tree);
use List::Util qw( min max );
use genomics::genotype;

# -----------------------
# -- Global variables  --
# -----------------------

our @ISA = qw(Exporter);
our @EXPORT = qw();

sub new
# Creation of new essay object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	my %tree;
	my %samples; # Samples info
	my %tmp_files; # TMP files to remove
	my %target_reference;
	my %modules;
	my @samples_order;
	my %references;
	
	my %defaults = (
		name => '',
		path => "/",
		queue => "shudra.q",
		priority => 0,
		pipeline => "",
		seq_platform => "",
		references => \%references,
		target_reference => \%target_reference,
		modules => \%modules,
		samples => \%samples,
		samples_order => \@samples_order,
		tree => \%tree,
		flow => {},
		flow_list => {},
		create_only_mode => 0,
		project_name => '',
	);
	my $self = {%defaults, %args};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub get_reads{
	my ($self,$sample,$replicate,$lane,$reads) = @_;
	return $self->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{$reads};
}

sub get_lane_type{
	my ($self,$sample,$replicate,$lane,$reads) = @_;
	return $self->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{type};
}

sub get_hold_jobs{
	my ($self,$module) = @_;
	my @previous_modules = split(",",$self->{flow}{$module}{previous});
	my $hold_jid = "";
	foreach(@previous_modules){
		my $state = $self->{modules}{$_}{state};
		if (substr($state,0,7) ne "pending"){
			my @tmp = split("=",$state);
			if (defined($tmp[1])){
				$hold_jid = $hold_jid . $tmp[1] . ",";
			}
		}
	}
	$hold_jid =~ s/,$//;
	return $hold_jid;
}

sub GetModuleHoldJobIDs{
	my ($self,$module) = @_;
	my @previous_modules = split(",",$self->{flow}{$module}{previous});
	my @jobs;
	foreach my $current_module (@previous_modules){
		foreach my $job (keys(%{$self->{modules}{$current_module}{jobs}})){
			push (@jobs,$self->GetJobID($current_module,$job));
		}
	}
	return @jobs;
}

sub GetModuleHoldString{
	my ($self,$module) = @_;
	my @jobs = $self->GetModuleHoldJobIDs($module);
	my $string = "";
	foreach my $job (@jobs){
		$string = $string . $self->GetSGEIDfromJobID($job) . ",";
	}
	if ($string ne ""){
		$string =~ s/,$//;
		$string = "-hold_jid ".$string
	}
	return $string;
}

sub GetSGEIDfromJobID{
	my ($self,$id) = @_;
	foreach my $module (keys(%{$self->{modules}})){
		foreach my $job (keys(%{$self->{modules}{$module}{jobs}})){
			if ($self->{modules}{$module}{jobs}{$job}{jobid} eq $id){
				return $self->{modules}{$module}{jobs}{$job}{SGE_id};
			}
		}
	}
	return "";
}

sub SetJobID{
	# Gets a new job ID, by looking all existent jobs numbers and assigning next number
	my ($self,$mymodule,$myjob) = @_;
	
	if (not defined($self->{modules}{$mymodule}{jobs}{$myjob}{jobid})){
		# Check all jobs in the different modules and retrieve last job number
		my $id = 0;
		foreach my $module (keys(%{$self->{modules}})){
			foreach my $job (keys(%{$self->{modules}{$module}{jobs}})){
				if (defined($self->{modules}{$module}{jobs}{$job}{jobid})){
					$id = max($id,$self->{modules}{$module}{jobs}{$job}{jobid});
				}
			}
		}
		$self->{modules}{$mymodule}{jobs}{$myjob}{jobid} = $id + 1;
	}
	return $self->{modules}{$mymodule}{jobs}{$myjob}{jobid};
}

sub GetJobID{
	# Returns job ID
	my ($self,$module,$job) = @_;
	return $self->{modules}{$module}{jobs}{$job}{jobid};
}

sub GetJobSGEID{
	my ($self,$module,$job) = @_;
	
	my $state = $self->GetJobState($module,$job);
	my $id = "";
	my @tmp = split("=",$state);
	if (defined($tmp[1])){
		$id = $tmp[1];
	} else {
		$id = $self->{modules}{$module}{jobs}{$job}{SGE_id};
	}
	return $id;
}

sub Precedes{
	# Deprecated function, previously used for sorting modules
	my ($modules,$mod1,$mod2) = @_;
	
	my @tmp = split(",",$$modules{$mod2}{previous});
	my $precedes = 0;
	if ($#tmp >=0){
		foreach my $prev (@tmp){
			if ($prev eq $mod1){
				$precedes =1;
			} else {
				$precedes = Precedes($modules,$mod1,$prev);
			}
		}
	}
	return $precedes;
}

sub GetLevel{
	# Function for getting the sorted level of a module within the whole analysis pipeline
	my ($self,$module) = @_;
	my $level = 0;
	if ($self->{flow}{$module}{previous} eq ""){
		$level = 1;
	} else {
		my @previous = split(",",$self->{flow}{$module}{previous});
		foreach my $mod (@previous){
			$level = max($self->GetLevel($mod) + 1,$level);
		}
	}
	$self->{flow}{$module}{level} = $level;
	return $level;
}

sub create_flow_list{
	my ($self) = @_;
	
	my @flow_list;
	my @modules = keys(%{$self->{flow}});
	foreach my $module (@modules){
		$self->GetLevel($module);
		if ($#flow_list < 0){
			push(@flow_list,$module);
		} else {
			my $place = 0;
			for (my $position=0;$position<=$#flow_list;$position ++){
				if ($self->{flow}{$flow_list[$position]}{level} < $self->{flow}{$module}{level}){
					$place = $position+1;
				}
			}
			splice(@flow_list,$place,0,$module);
		}
	}
	
	$self->{flow_list} = \@flow_list;
}

sub checks{
	my ($essay, $references) = @_;
	my $specie = "";
	my $seq_platform = "";
	my $capture_system = "";
	foreach my $sample (@{$essay->{samples_order}}){
		# Check specie
		if ($specie eq ""){
			$specie = $essay->{'samples'}{$sample}{'specie'};
		} elsif ($essay->{samples}{$sample}{specie} ne $specie){
			die "ERROR: Specie $essay->{samples}{$sample}{specie} for $sample not equal to essay specie $specie!! Please check species!!\n";
		}
		
		# Check sequencing platform & capture system
		foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
			my $current_seq_platform = $essay->{samples}{$sample}{replicates}{$replicate}{'seq_platform'};
			my $current_capture_system = $essay->{samples}{$sample}{replicates}{$replicate}{'capture_system'};
			if ($seq_platform eq ""){
				$seq_platform = $current_seq_platform;
			} elsif ($current_seq_platform ne $seq_platform){
				$seq_platform = "mixed";
			}
			if ($capture_system eq ""){
				$capture_system = $current_capture_system;
			} elsif ($current_capture_system ne $capture_system){
				$capture_system = "mixed";
			}
			# Check and store read lengths
			foreach my $lane (sort(keys(%{$essay->{samples}{$sample}{replicates}{$replicate}{lanes}}))){
				if ($current_seq_platform eq "Illumina"){
					my $R1 = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{R1};
					my $R2 = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{R2};
					$essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{length_R1} = `gunzip -c $R1 | head -2 | tail -1 | wc | awk '{print \$3}'` - 2;
					$essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{length_R2} = `gunzip -c $R2 | head -2 | tail -1 | wc | awk '{print \$3}'` - 2;
				} elsif ($current_seq_platform eq "SOLiD"){
					my $F3 = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F3_csfasta};
					my $F5 = $essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{F5_csfasta};
					$essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{length_R1} = `head -2 $F3 | tail -1 | wc | awk '{print \$3}'` - 2;
					$essay->{samples}{$sample}{replicates}{$replicate}{lanes}{$lane}{data}{length_R1} = `head -2 $F5 | tail -1 | wc | awk '{print \$3}'` - 2;
				} else {
					die "ERROR: Invalid sequencing platform $current_seq_platform!!\n";
				}
			}
		}
		
		# Check capture system
		
	}
	$essay->{seq_platform} = $seq_platform;
	
	# Check essay type (pipeline). Be sure to add here any pipeline specific conditions for retrieving information from config json:
	my $tmp;
	if (substr($essay->{pipeline},0,9) eq "DNA-reseq" or (substr($essay->{'pipeline'},0,2) eq "WG")){
		$tmp = $$references{genomes};
		my $found = 0;
		foreach (keys(%{$tmp})){
			if ($specie eq $_){
				$found = 1;
			}
		}
		if ($found){
			$essay->{references} = $$references{genomes}{$specie};
		} else {
			die "ERROR: Samples specie ($specie) isn't included in project references!!\n";
		}
	}
	
	# Load pipeline information. Be sure that pipelines are executable and can be found in path
	my $pipeline_json = which($essay->{pipeline}.".json");
	if ($pipeline_json eq ""){
		die "ERROR: Pipeline $essay->{pipeline} not available!!\n";
	}
	my %pipe = get_from_json($pipeline_json);
	
	# Store pipeline flow into essay object
	$essay->{flow} = $pipe{modules};
	
	# Create and store in essay the list of modules (sorted to be launched)
	$essay->create_flow_list();
}

sub is_solid{
	my ($essay) = @_;
	if ($essay->{seq_platform} eq "SOLiD"){
		return 'true';
	} else {
		return 'false';
	}
}

sub exists_module{
	my ($essay, $module) = @_;
	if (defined $essay->{flow}{$module}){
		return 1;
	} else {
		return 0;
	}
}

sub get_samples {
	my ($self, $samples) = @_;
	my $n = 1;
	foreach my $sample (&get_keys (\$samples)) {
		foreach my $replicate (&get_keys (\$$samples{$sample}{replicates}{ids})) {
			foreach my $lane (&get_keys (\$$samples{$sample}{replicates}{ids}{$replicate}{lanes}{ids})) {
				$self->{samples}{$n} = lane->get_lane ($$samples{$sample}{replicates}{ids}{$replicate}{lanes}{ids}{$lane});
				$self->{files2merge} .= ",$self->{samples}{$n}->{rg_file}";
				$self->{readgroups} .= "$self->{samples}{$n}->{rgsm},";
				$n++;
			}
		}
	}
	$self->{readgroups} =~ s/,$//;
	return 1;
}

sub create_dirs_to_delete {
	my ($tree, $path) = @_;
	foreach my $folder (&get_keys($tree)){
		if (ref($$tree{$folder}) eq "HASH" or $$tree{$folder} eq ''){
			# First, new folder is created
			make_path ($path."/".$folder, {verbose => 1, mode => 0755,});
		}
		if (ref($$tree{$folder}) eq "HASH"){
			# Subfolders of current folder must be also created
			create_dirs($$tree{$folder},$path."/".$folder)
		}
	}	
}

sub create_dirs {
	my ($tree, $path) = @_;
	foreach my $folder (&get_keys($tree)){
		if (ref($$tree{$folder}) eq "HASH"){
			# First, new folder is created
			make_path ($path."/".$folder, {verbose => 1, mode => 0755,});
			# Subfolders of current folder must be also created
			create_dirs($$tree{$folder},$path."/".$folder)
		}
	}	
}

sub create_tree {
	my ($self) = @_;
	# if essay root path does not exist, create it
	unless(-d $self->{path}){
		make_path ($self->{path}, {verbose => 1, mode => 0755,});
	}
	&create_dirs ($self->{tree},$self->{path});
	return 1;
}

sub add_element_to_hash{
	my ($self,$name,$value,@keys) = @_;
	if ($#keys > 0){
		my $key = shift(@keys);
		$$self{$key} = add_element_to_hash($$self{$key},$name,$value,@keys);
	} elsif ($#keys == 0) {
		if ($name eq ""){
			$$self{$keys[0]} = {};
		} else {
			$$self{$keys[0]}{$name} = $value;
		}
	} else {
		$$self{$name} = $value;
	}
	return $self;
}

sub add_file{
	my ($self,@params) = @_;
	my $path = $self->create_path(@params);
	my $name = pop(@params);
	add_element_to_hash($self->{tree},$name,$path,@params);
	return $self->get_path(@params,$name);
}

sub mkdir{
	my ($self,@folders) = @_;
	add_element_to_hash($self->{tree},"","",@folders);
	return $self->get_path(@folders);
}

sub get_path{
	my ($self,@folders) = @_;
	
	# Check if file or directory is already in the tree structure
	if (not defined(get_hash_parameter($self->{tree},@folders))){
		return undef;
	} else {
		return $self->create_path(@folders);
	}
}

sub create_path{
	my ($self,@folders) = @_;
	my $path = "";
	if ($#folders > 0){
		my $name = pop(@folders);
		$path = $self->create_path(@folders) . "/" . $name;
	} else {
		$path = $self->{path} . "/" . $folders[0];
	}
	# Check if tag is already present in hash
	if (ref(get_hash_parameter($self->{tree},@folders)) eq ""){
		if (defined(get_hash_parameter($self->{tree},@folders))){
			$path = get_hash_parameter($self->{tree},@folders);
		}
	}
	return $path;
}

sub get_hash_parameter{
	my ($self,@keys) = @_;
	my $result;
	if ($#keys > 0){
		my $key = shift(@keys);
		$result = &get_hash_parameter($$self{$key},@keys);
	} else {
		$result = $$self{$keys[0]};
	}
	return $result;
}

sub get_parameter{
	my ($self,@keys) = @_;
	my $result;
	if ($#keys > 0){
		my $key = shift(@keys);
		$result = &get_hash_parameter($self->{$key},@keys);
	} else {
		$result = $self->{$keys[0]};
	}
	return $result;
}

sub get_readgroups {
	my ($essay,$type) = @_;
	my @readgroups;
	foreach my $sample (@{$essay->{samples_order}}){
		foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
			foreach my $lane (sort(keys(%{$essay->{samples}{$sample}{replicates}{$replicate}{lanes}}))){
				if ($essay->is_solid() eq 'true'){
					if (defined($type)){
						push (@readgroups,$sample.".".$replicate.".".$lane.".".$type);
					} else {
						die "ERROR: No type (PE or F3) defined por SOLiD analysis!!\n";
					}
				} else {
					push (@readgroups,$sample.".".$replicate.".".$lane);
				}				
			}
		}
	}
	return \@readgroups;
}

sub get_replicates {
	my ($essay,$type) = @_;
	my @replicates;
	foreach my $sample (@{$essay->{samples_order}}){
		foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
			if (defined($type)){
				push (@replicates,$sample.".".$replicate.".".$type);
			} else {
				push (@replicates,$sample.".".$replicate);
			}
		}
	}
	return \@replicates;
}

sub get_bams {
	# This function returns an array containing the bams of a sample, replicate, lane, etc.
	# Consider that there are bams for each lane and replicate, even chromosome-splitted, 
	# and at different stages of the analysis (Q1-filtering, realignment, etc).
	# - Use parameter "level" for defining which bam-level you desire: "sample" (if existent),
	#   "replicate" or "lane" level
	# - Use parameter "chr" for requesting chromosome-splitted bams: "all" for all chromosomes,
	#   or a value in (1..22,X,Y,M) for a given chromosome-splitted bam
	# - Use parameter "stage" for requesting bam at a given stage (not available at the moment)
	my $essay = shift;
	my %defaults = (sample=>undef,replicate=>undef,lane=>undef,type=>undef,level=>"lane",chr=>undef);
	my $parameters = {%defaults,@_};
	
	my @lanes_bams;
	foreach my $sample (@{$essay->{samples_order}}){
		if ((not defined($$parameters{sample})) or ($sample eq $$parameters{sample})) {
			if ($$parameters{level} ne "sample"){
				foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
					if ((not defined($$parameters{replicate})) or ($replicate eq $$parameters{replicate})){
						if ($$parameters{level} ne "replicate"){
							foreach my $lane (sort(keys(%{$essay->{samples}{$sample}{replicates}{$replicate}{lanes}}))){
								if ((not defined($$parameters{lane})) or ($lane eq $$parameters{lane})) {
									if ($essay->is_solid() eq 'true'){
										if ($$parameters{type} eq "PE"){
											push(@lanes_bams,$essay->get_path("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".readgroup.bam"));
										} else {
											push(@lanes_bams,$essay->get_path("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane."_".$$parameters{type}.".readgroup.bam"));
										}
									} else {
										push(@lanes_bams,$essay->get_path("trash","mapping",$sample,$replicate,$lane,$sample."_".$replicate."_".$lane.".readgroup.bam"));
									}
								}
							}
						} else {
							if (defined($$parameters{chr})){
								foreach my $chr (@{genotype->Chromosomes()}){
									if ($chr eq "all" or $chr eq $$parameters{chr}){
										if ($essay->is_solid() eq 'true' and $$parameters{type} eq "F3"){
											push(@lanes_bams,$essay->get_path("trash","variants",$sample,$replicate,$sample."_".$replicate."_split".$$parameters{type}).".REF_chr".$chr.".bam");
										} else {
											push(@lanes_bams,$essay->get_path("trash","variants",$sample,$replicate,$sample."_".$replicate."_split").".REF_chr".$chr.".bam");
										}
									}
								}
							} else {
								if ($essay->is_solid() eq 'true'){
									if ($$parameters{type} eq "PE"){
										push(@lanes_bams,$essay->get_path("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam"));
									} else {
										push(@lanes_bams,$essay->get_path("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_Q1_PCR_realign_F3.bam"));
									}
								} else {
									push(@lanes_bams,$essay->get_path("analysis","mapping",$sample,$replicate,$sample."_".$replicate.".bam"));
								}
							}
						}
					}
				}
			} else {
				die "ERROR: No sample-level bam has being created!\n";
			}
		}
	}
	return \@lanes_bams;
}

sub get_keys {
	my ($self) = @_;
	return (sort {$a cmp $b} keys %{$self});
}

sub GetJobState {
	my ($self,$module,$job) = @_;
	if (defined($self->{modules}{$module}{jobs}{$job}{state})){
		return $self->{modules}{$module}{jobs}{$job}{state};
	} else {
		return "pending";
	}
}

sub GetModuleState {
	my ($self,$module) = @_;
	if (defined($self->{modules}{$module}{state})){
		return $self->{modules}{$module}{state};
	} else {
		return "";
	}
}

sub SetJobState {
	my ($self,$module,$job,$state) = @_;
	$self->{modules}{$module}{jobs}{$job}{state} = $state;
}

sub SetJobSGEID {
	my ($self,$module,$job,$sgeid) = @_;
	$self->{modules}{$module}{jobs}{$job}{SGE_id} = $sgeid;
}

sub SetJobLog {
	my ($self,$module,$job,$log) = @_;
	$self->{modules}{$module}{jobs}{$job}{log} = $log;
}

sub SetHoldJobs {
	my ($self,$module,$job,$hold) = @_;
	$self->{modules}{$module}{jobs}{$job}{hold_jobs} = $hold;
}

sub update {
	my ($self) = @_;
	create_json($self,$self->{path} . "/state");
}

sub get_state {
	my ($self) = @_;
	my %tmp = &get_from_json($self->{path} . "/state.json");
	foreach my $key (keys(%tmp)){
		$self->{$key} = $tmp{$key};
	}
	$self = %tmp;
}

1;