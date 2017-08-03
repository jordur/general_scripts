#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with bioinformatics analysis information
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;

# -----------------------
# -- Global variables  --
# -----------------------

package variants_info;
sub new
# Creation of new variants_information object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my %info = (
		"filter_col" => 7,
		"samples" => undef, # It's very important to start with the "undef" value, otherwise problems arise!!!!
		"Strand_Bias" => 0,
	);
	
	my %format_data = (
		"format_col" => 8,
		"data_col" => 9,
		"ids" => undef, # It's very important to start with the "undef" value, otherwise problems arise!!!!
	);
	my %defaults = (
		Chr => 0,
		Pos => 1,
		ID => 2,
		Ref_Allele => 3,
		Var_Allele => 4,
		qual => 5,
		filter => 6,
		info => \%info,
		format_data => \%format_data,
	);
	my $self = {%defaults,@_};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub set_parameter ()
# Function for setting parameters into an already created analysis_info object
{
	my ($self, $param, $value) = @_;
	if ($param eq "gatk" or $param eq "GK"){
		$self->{format_data}->{ids}->{GK} = $value;
	} elsif ($param eq "varscan" or $param eq "VS"){
		$self->{format_data}->{ids}->{VS} = $value;
	} elsif ($param eq "PE" or $param eq "PA"){
		$self->{format_data}->{ids}->{PA} = $value;
	} elsif ($param eq "F3"){
		$self->{format_data}->{ids}->{F3} = $value;
	} elsif ($param eq "GFF" or $param eq "GF"){
		$self->{format_data}->{ids}->{GF} = $value;
	} elsif ($param eq "GT"){
		$self->{format_data}->{ids}->{GT} = $value;
	} elsif (defined($self->{$param})){
		$self->{$param} = $value;
	} elsif (defined($self->{format_data}->{$param})){
		$self->{format_data}->{$param} = $value;
	} elsif (defined($self->{info}->{$param})){
		$self->{info}->{$param} = $value;
	}
}

sub get_parameter ()
# Function for getting parameters from analysis_info objects
{
	my ($self, $param) = @_;
	if ($param eq "gatk" or $param eq "GK"){
		return $self->{format_data}->{ids}->{GK};
	} elsif ($param eq "varscan" or $param eq "VS"){
		return $self->{format_data}->{ids}->{VS};
	} elsif ($param eq "PE" or $param eq "PA"){
		return $self->{format_data}->{ids}->{PA};
	} elsif ($param eq "F3"){
		return $self->{format_data}->{ids}->{F3};
	} elsif ($param eq "GFF" or $param eq "GF"){
		return $self->{format_data}->{ids}->{GF};
	} elsif ($param eq "GT"){
		return $self->{format_data}->{ids}->{GT};		
	} elsif (defined($self->{$param})){
		return $self->{$param};
	} elsif (defined($self->{info}->{$param})){
		return $self->{info}->{$param};
	} elsif (defined($self->{format_data}->{$param})){
		return $self->{format_data}->{$param};
	} else {
		return undef;
	}
}

sub set_data_parameter ()
# Function for setting parameters from data section of variants_info objects
{
	my ($self, $param, $value) = @_;
	$self->{format_data}->{ids}->{$param} = $value;
}

sub get_data_parameter ()
# Function for getting parameters from data section of variants_info objects
{
	my ($self, $param) = @_;
	return $self->{format_data}->{ids}->{$param};
}

sub set_sample_parameter ()
# Function for setting parameters from samples of variants_info objects
{
	my ($self, $sample, $param, $value) = @_;
	$self->{info}->{samples}->{$sample}->{$param} = $value;
}

sub get_sample_parameter ()
# Function for getting parameters from samples of variants_info objects
{
	my ($self, $sample, $param) = @_;
	if (defined($self->{info}->{samples}->{$sample}->{$param})){
		return $self->{info}->{samples}->{$sample}->{$param};
	} else {
		return;
	}
}

sub get_ids ()
{
	my ($self) = @_;
	return keys (%{$self});
}

sub get_data_ids ()
{
	my ($self) = @_;
	return keys (%{$self->{format_data}->{ids}});
}

sub get_sample_ids ()
{
	my ($self) = @_;
	return keys(%{$self->{info}->{samples}});
}

sub get_fields_from_vcf
# Function that retrieves information fields from vcf files
{
	my ($class, $vcf_line, $config) = @_;
	my $self = variants_info->new();
	
	# Retrieve information
	my @vcf_cols = split ('\t', $vcf_line);
	my @variant_ids = $config->get_ids();
	foreach (@variant_ids){
		my $config_key = $_;
		if ($config_key ne "info" and $config_key ne "format_data"){
			$self->set_parameter($config_key,$vcf_cols[$config->{$config_key}]);
		}
	}
		
	# Information
	my @cols = split (';', $vcf_cols[$config->get_parameter("filter_col")]);
	my %samples;
	my @sample_ids = $config->get_sample_ids();
	foreach (@sample_ids){
		my $sample_name = $_;
		my @tmp = split('=', $cols[$config->get_sample_parameter($sample_name,"Depth")]);
		$self->set_sample_parameter($sample_name,"Depth",$tmp[1]);
		@tmp = split('=', $cols[$config->get_sample_parameter($sample_name,"Freq")]);
		$self->set_sample_parameter($sample_name,"Freq",$tmp[1]);
		@tmp = split('=', $cols[$config->get_sample_parameter($sample_name,"Genotype")]);
		$self->set_sample_parameter($sample_name,"Genotype",$tmp[1]);
	}
	my @tmp = split ('=', $cols[$config->get_parameter("Strand_Bias")]);
	$self->set_parameter("Strand_Bias",$tmp[1]);
	
	# format_data information
	@cols = split (':', $vcf_cols[$config->get_parameter("data_col")]);
	my @data_ids = $config->get_data_ids();
	foreach (@data_ids){
		my $id_name = $_;
		$self->set_data_parameter($id_name,$cols[$config->get_parameter($id_name)]);
	}
		
	# Return object
	return $self;
}

sub put_fields_into_vcf
# Function that returns a vcf line containing the information in an variants_info object
{
	my ($data, $config) = @_;
	my $line;
	$line = $data->get_parameter("Chr") . "\t" .  $data->get_parameter("Pos") . "\t" . $data->get_parameter("ID") . "\t" . $data->get_parameter("Ref_Allele") . "\t" .
		$data->get_parameter("Var_Allele") . "\t" . $data->get_parameter("qual") . "\t" . $data->get_parameter("filter") . "\t";

	# Filter info:
	my @samples_ids = $config->get_sample_ids();
	for (my $i=1;$i<=$#samples_ids+1;$i++){
		my $sample_key = "S" . $i;
		if ($i > 1){
			$line .= ";";
		}
		$line .= $sample_key . "D=" . $data->get_sample_parameter($sample_key,"Depth") . ";" .
			$sample_key . "F=" . $data->get_sample_parameter($sample_key,"Freq") . ";" . 
			$sample_key . "G=" . $data->get_sample_parameter($sample_key,"Genotype");
	}
	$line .= ";SB=".$data->get_parameter("Strand_Bias")."\t";
	
	# Format_data:
	my $format = "";
	my $data2 ="";
	my @data_ids = $config->get_data_ids();
	for (my $i=0;$i<=$#data_ids;$i++){
		my $id_key = "";
		foreach (@data_ids){
			my $data_id = $_;
			if ($config->get_parameter($data_id) == $i){
				$id_key = $data_id;
			}
		}
		if ($i > 0){
			$format .= ":";
			$data2 .= ":";
		}
		$format .= $id_key;
		$data2 .= $data->get_parameter($id_key);
	}
	$line .= $format . "\t" . $data2 . "\n";
}

sub create_from_vcf
# Creation of new variants_information object from a vcf file
{
	my ($class, $vcf_file) = @_;
	open (VCF, "<", $vcf_file) or die("ERROR: Couldn't open file $vcf_file!! $!\n");
	my $self = variants_info->new();
	
	# Get header (first line of file) and check the amount of samples in the file:
	my $header = <VCF>;
	my @vcf_cols = split ('\t', $header);
	my %info;
	my %format_data;
	$info{filter_col}=7;
	$format_data{format_col}=8;
	$format_data{data_col}=9;
	
	# Third line of file is needed for checking samples number, data format, etc...
	my $second_line = <VCF>;
	my $third_line = <VCF>;
	close VCF;
	
	if (defined ($third_line)){
		@vcf_cols = split ('\t', $third_line);
	
		# information
		my @cols = split (';', $vcf_cols[$info{filter_col}]);
		my $col_number = 0;
		my $sample_num = 0;
		foreach (my $i=0;$i<=$#cols;$i++){
			my $col = $cols[$i];
			my @tmp = split ("=", $col);
			if ($tmp[0] ne "SB"){
				my $sample_name = "S" . ($sample_num + 1);
				$self->set_sample_parameter($sample_name,"Depth",$i);
				$self->set_sample_parameter($sample_name,"Freq",$i + 1);
				$self->set_sample_parameter($sample_name,"Genotype",$i + 2);
				$i = $i + 2;
				$sample_num ++;
			} else {
				$self->set_parameter("Strand_Bias",$i);
			}
		}
		
		# format_data information
		@cols = split (':', $vcf_cols[$format_data{format_col}]);
		$col_number = 0;
		foreach (@cols){
			my $col = $_;
			$self->set_parameter($col,$col_number);
			$col_number++;
		}
	} else {
		die "ERROR: Incorrect or empty file $vcf_file!!!\n";
	}
		
    # Return object
    return $self;
}

sub get_config () 
# Get variants_information configuration of a variants vcf
{
	my $class = shift;
	my %defaults = (
		vcf_file => '',
	);
	my $attributes = {%defaults, @_};
	
	my $values = variants_info->create_from_vcf($attributes->{vcf_file});
	return $values;
}

1;
