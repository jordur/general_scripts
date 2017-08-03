#!/usr/bin/perl -w
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with bioinformatics analysis information
# @Author: Arbol & JM Rosa
# @Contributors:
########################################################

package analysis_info;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use Scalar::Util;
use Storable qw(dclone);
use genomics::genotype;

# -----------------------
# -- Global variables  --
# -----------------------
my $columns_number = 38;

sub new
# Creation of new analysis_information object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my %defaults = (
		HGNC_symbol => 0,
		Chr => 1,
		Pos => 2,
		Ref_Allele => 3,
		Var_Allele => 4,
		Strand_Bias => '',
		Gene => '',
		Gene_description => '',
		HGVSc_name => '',
		Intron => '',
		Exon => '',
		HGVSp_name => '',
		Variant_effect => '',
		ALL_MAF => '', 
		AFR_MAF => '',
		AMR_MAF => '',
		ASN_MAF => '',
		EUR_MAF => '',
		InterPro_IDs => '',
		InterPro_descriptions => '',
		HGMD_info => '',
		Related_publication => '',
		Existing_variation => '',
		Feature_type => '',
		Feature_ID => '',
		RefSeq_ID => '',
		CCDS_ID => '',
		Canonical_isoform => '',
		Conservation_score => '',
		Grantham_distance => '',
		Condel_prediction => '',
		SIFT_prediction => '',
		PolyPhen_prediction => '',
		Affected_prot_domains => '',
		Regulatory_Motif_name => '',
		Regulatory_Motif_position => '',
		Regulatory_High_Inf_Pos => '',
		Regulatory_Motif_Score_Change => '',
		Flanking_sequence => '',
		type => 'both',
		separator => '|',
		samples => '',
		print_known => 'true',
		print_unknown => 'true',
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
	$self->{$param} = $value;
}

sub get_parameter ()
# Function for getting parameters from analysis_info objects
{
	my ($self, $param) = @_;
	if (defined($self->{$param})){
		return $self->{$param};
	} else {
		print Dumper($self);
		die "ERROR: No parameter $param was found in analysis_info object!!\n";
	}
}

sub set_sample_parameter ()
# Function for setting parameters from samples of variants_info objects
{
	my ($self, $sample, $param, $value) = @_;
	$self->{samples}->{$sample}->{$param} = $value;
}

sub get_sample_parameter ()
# Function for getting parameters from samples of analysis_info objects
{
	my ($self, $sample, $param) = @_;
	if (defined($self->{samples}->{$sample}->{$param})){
		return $self->{samples}->{$sample}->{$param};
	} else {
		return;
	}
}

sub get_sample_ids ()
{
	my ($self) = @_;
	return keys(%{$self->{samples}});
}

sub compare()
# Function to compare an analysis_info object to another. Chr, position, ref_allele and var_allele are checked. If equal, 0 is returned. If not, 1 is returned if $self is "before" $comp_var, -1 if the other case.
{
	my ($self,$comp_var) = @_;
	if (genotype->GetChromosomeNumber($self->get_parameter("Chr")) > genotype->GetChromosomeNumber($comp_var->get_parameter("Chr"))){
		return -1;
	} elsif (genotype->GetChromosomeNumber($self->get_parameter("Chr")) < genotype->GetChromosomeNumber($comp_var->get_parameter("Chr"))){
		return 1;
	} elsif (int($self->get_parameter("Pos")) > int($comp_var->get_parameter("Pos"))){
		return -1;
	} elsif (int($self->get_parameter("Pos")) < int($comp_var->get_parameter("Pos"))){
		return 1;
	} elsif ($self->get_parameter("Ref_Allele") gt $comp_var->get_parameter("Ref_Allele")){
		return -1;
	} elsif ($self->get_parameter("Ref_Allele") lt $comp_var->get_parameter("Ref_Allele")){
		return 1;
	} elsif ($self->get_parameter("Var_Allele") gt $comp_var->get_parameter("Var_Allele")){
		return -1;
	} elsif ($self->get_parameter("Var_Allele") lt $comp_var->get_parameter("Var_Allele")){
		return 1;
	} else {
		return 0;
	}
}

sub get_fields_from_psv ()
# Function that returns the values of a line of a psv file into a analysis_info objetct
{
	my ($class, $line, $config) = @_;
	my @data = split (/[$config->{separator}]/, $line);
	my $temp = analysis_info->new();
	
	foreach my $key_name (keys (%{$config})){
		if($key_name eq "samples"){
			my %samp;
			foreach (keys (%{$config->{$key_name}})){
				my $sample = $_;
				$samp{$sample}{Depth} = $data[$config->{$key_name}->{$sample}->{Depth}];
				$samp{$sample}{Freq} = $data[$config->{$key_name}->{$sample}->{Freq}];
				$samp{$sample}{Genotype} = $data[$config->{$key_name}->{$sample}->{Genotype}];
			}
			$temp->set_parameter("samples",\%samp);
		} else {
			if (Scalar::Util::looks_like_number($config->{$key_name})){
				$temp->set_parameter($key_name, $data[$config->{$key_name}]);
			}
		}
	}
	# Return object
	return $temp;
}

sub put_fields_into_psv ()
# Function that returns a tsv line containing the fields of an analysis info objetct
{
	my ($data, $config) = @_;
	my $line = "";
	my $sample_num=scalar(keys(%{$config->{samples}}));
	for(my $i=0;$i<=($sample_num*3)+$columns_number;$i++){
		if ($i<5 or $i>($sample_num*3)+4){
			foreach my $key (keys(%{$config})){
				if ($config->{$key} eq $i){
					$line .= $data->{$key};
					next;
				}
			}
		} else {
			foreach my $sample (@{$config->{samples_order}}){
				foreach my $sample_field (keys(%{$config->{samples}->{$sample}})){
					if ($config->{samples}->{$sample}->{$sample_field} eq $i){
						$line .= $data->{samples}->{$sample}->{$sample_field};
						next;
					}
				}
			}
		}
		$line .= $config->{separator};
	}
	chop($line);
	$line.="\n";
	return $line;
}

sub create_from_psv
# Creation of new analysis_information object from a psv file
{
	my ($class, $type, $psv_file, $sep_id) = @_;
	open (PSV, "<", $psv_file) or die("ERROR: Couldn't open file $psv_file!! $!\n");
	
	# Get header (first line of file) and check the amount of samples in the file:
	my $header = <PSV>;
	my @cols = split (/[$sep_id]/, $header);
	my %samples;
	my @samples_order;
	foreach (@cols){
		my $col = $_;
		my $tmp;
		if (length($col) > 9){
			$tmp = substr($col, -9);
		} else {
			$tmp = $col;
		}
		if ($tmp eq "_Genotype"){
			my $sample_name = substr($col,0,length($col) - 9);
			push (@samples_order,$sample_name);
		}
	}
	
	my $self = analysis_info->create_config(\@samples_order,$type,$sep_id);    
    close PSV;
    
    # Return object
    return $self;
}

sub GetGeneName
# Function that returns a suitable gene name (either HGNC, or ENSEMBL ID, ...)
{
	my ($self) = @_;
	# Define gene name (not all genes have HGNC_symbol). If no information available, report variant as regulatory_region
	my $gene; # Remember that "regulatory_region" is the reserved name for those variants in UTRs or regulatory regions
	if ($self->{HGNC_symbol} ne '-'){
		$gene = $self->{HGNC_symbol};
	} elsif ($self->{Gene} ne '-'){
		$gene = $self->{Gene};
	} elsif ($self->{Feature_ID} ne '-'){
		$gene = $self->{Feature_ID};
	} else {
		$gene = "unknown_regulatory_region";
	}
	return $gene;
}

sub create_config
# Create config object from samples, type and separator
{
	my ($class,$samples_list,$type,$sep_id) = @_;
	my $sample_num = $#{$samples_list} + 1;
	my %samples;

	# Define values in samples hash
	for (my $i=0;$i<=$#{$samples_list};$i++){
		$samples{$$samples_list[$i]}{Genotype} = $i*3+5;
		$samples{$$samples_list[$i]}{Depth} = $i*3+6;
		$samples{$$samples_list[$i]}{Freq} = $i*3+7;
	}
	
	my $new_config = analysis_info->new(
		HGNC_symbol => 0,
		Chr => 1,
		Pos => 2,
		Ref_Allele => 3,
		Var_Allele => 4,
		Strand_Bias => ($sample_num*3) + 5,
		Gene => ($sample_num*3) + 6,
		Gene_description => ($sample_num*3) + 7,
		HGVSc_name => ($sample_num*3) + 8,
		Intron => ($sample_num*3) + 9,
		Exon => ($sample_num*3) + 10,
		HGVSp_name => ($sample_num*3) + 11,
		Variant_effect => ($sample_num*3) + 12,
		ALL_MAF => ($sample_num*3) + 13, 
		AFR_MAF => ($sample_num*3) + 14,
		AMR_MAF => ($sample_num*3) + 15,
		ASN_MAF => ($sample_num*3) + 16,
		EUR_MAF => ($sample_num*3) + 17,
		InterPro_IDs => ($sample_num*3) + 18,
		InterPro_descriptions => ($sample_num*3) + 19,
		HGMD_info => ($sample_num*3) + 20,
		Related_publication => ($sample_num*3) + 21,
		Existing_variation => ($sample_num*3) + 22,
		Feature_type => ($sample_num*3) + 23,
		Feature_ID => ($sample_num*3) + 24,
		RefSeq_ID => ($sample_num*3) + 25,
		CCDS_ID => ($sample_num*3) + 26,
		Canonical_isoform => ($sample_num*3) + 27,
		Conservation_score => ($sample_num*3) + 28,
		Grantham_distance => ($sample_num*3) + 29,
		Condel_prediction => ($sample_num*3) + 30,
		SIFT_prediction => ($sample_num*3) + 31,
		PolyPhen_prediction => ($sample_num*3) + 32,
		Affected_prot_domains => ($sample_num*3) + 33,
		Regulatory_Motif_name => ($sample_num*3) + 34,
		Regulatory_Motif_position => ($sample_num*3) + 35,
		Regulatory_High_Inf_Pos => ($sample_num*3) + 36,
		Regulatory_Motif_Score_Change => ($sample_num*3) + 37,
		Flanking_sequence => ($sample_num*3) + 38,
		type => $type,
		separator => $sep_id,
		samples => \%samples,
		samples_order => $samples_list,
    );
}

sub get_config () 
# Get analysis_information configuration of a results psv
{
	my $class = shift;
	my %defaults = (
		type => 'both',
		sample_num => 1,
		sep_id => '|',
		psv_file => '',
	);
	my $attributes = {%defaults, @_};
	
	my $values = analysis_info->create_from_psv($attributes->{type}, $attributes->{psv_file}, $attributes->{sep_id});
	if ($attributes->{type} eq 'known'){
		$values -> {type} = 'known';
		$values -> {print_unknown} = 'false';
	}
	if ($attributes->{type} eq 'unknown'){
		$values -> {type} = 'unknown';
		$values -> {print_known} = 'false';
	}
	return $values;
}

sub create_header () 
# Create psv file header from configuration object
{
	my $config = shift;
	
	my $line = "#";
	my $sample_num=scalar(keys(%{$config->{samples}}));
	for(my $i=0;$i<=($sample_num*3)+$columns_number;$i++){
		if ($i<5 or $i>($sample_num*3)+4){
			foreach my $key (keys(%{$config})){
				if ($config->{$key} eq $i){
					$line .= $key;
					next;
				}
			}
		} else {
			# TODO: In the future, "_Var/Depth" should be also changed to "Freq"!!!
			foreach my $sample (@{$config->{samples_order}}){
				foreach my $sample_field (keys(%{$config->{samples}->{$sample}})){
					if ($config->{samples}->{$sample}->{$sample_field} eq $i){
						if ($sample_field eq "Freq"){
							$line .= $sample."_Var/Depth";
						} else {
							$line .= $sample."_".$sample_field;
						}
						next;
					}
				}
			}
		}
		$line .= $config->{separator};
	}
	chop($line);
	return $line;
}

1;