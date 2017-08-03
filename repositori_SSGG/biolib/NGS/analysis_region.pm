#!/usr/bin/perl -w
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with CNVs
# @Author: JM Rosa
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


package analysis_region;
sub new
# Creation of new analysis_information object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my %defaults = (
		chrom => '',
		c_chr => '',
		start => '',
		end => '',
		value => '',
	);

	my $self = {%defaults,@_};
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub get_c_chr {
#Get a numerical value for the chromosome	
	my $chr = shift;
	$chr =~ s/chr//g;
	
	my $result;
	if ($chr eq 'X'){
		$result = 23;
	}
	elsif ($chr eq 'Y'){
		$result = 24;
	}
	elsif ($chr eq 'M'){
		$result = 25;
	}
	else {
		$result = $chr;
	}
	return $result;
}

sub get_region () 
# Get the region
{
	my $class = shift;
	my $info = shift;
	
	my @data = split (/\s+/, $info);
	
	my $c_chr = &get_c_chr ($data[3]);
	
	my $values = analysis_region->new (
		chrom => $data[3], 
		c_chr => $c_chr, 
		start => $data[4], 
		end => $data[5], 
		value => $data[7]
	);
	
	# Return object
	return $values;
}

sub get_log_ratio {
	my $clas = shift;
	my ($tumour, $normal) = @_;
	
	my $t_object = analysis_region->get_region($tumour);
	my $n_object = analysis_region->get_region($normal);
	
	my $value = &log10 ( (exp ($t_object->{value}) / (exp ($n_object->{value}) ) ) );
	#print "$value\n";
	return $value;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

1;

