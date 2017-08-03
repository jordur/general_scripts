#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Exclude Variants present in a list of samples
# @Author: Arbol
# @Contributors:
########################################################

use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;
use File::Copy;
use File::Path;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;
use analysis_info;
use file_handle;
use Scalar::Util;

my $input;
my $sep = '|';
my $output = "excluded.psv";
my $samples;
my @samples;

GetOptions (	
				"i=s"	=> \$input,	
				"d=s"	=> \$sep,
				"o=s"	=> \$output,	
				"s=s"	=> \$samples,	
);

print "\n\n=================================================================
Script for excluding variants present in a list of samples, from SG output annotated variants files :: v 0.1\n";

if (not defined $input and not defined $samples) {
	die "
			Parameters:
		-i file_name        Input file name (both SG anotated files -csv- or variants files -vcf- are supported)
		-d separator        Field separator (default = '|')
		-o output_name      Output name (optional, default name is \"excluded.psv\")
		-s samples          Comma separated list of sample names
\n".localtime()."\n=================================================================\n\n";
}

&main ();
print  "\nProcess finished... ".localtime()."\n";

sub main () {
	# Config will be loaded, ouput file handles are created
	my $config = analysis_info->get_config(psv_file => $input);
	my $output_handle = &get_out_file_handle($output);
	
	# Creation of samples array from samples csv
	@samples = split(",",$samples);
	
	my $input_handle = &get_in_file_handle($input);
	while (my $line = <$input_handle>){
		chomp $line;
		if ($line =~ m/^#/) {
			print $output_handle "$line\n";
		} else {
			if (not &exclude($line,$config)){
				print $output_handle "$line\n";
			}
		}
	}
}

sub IsVariant(){
	my ($genotype) = @_;
	if ($genotype eq "UNC_Hetero" or $genotype eq "P_Hetero" or $genotype eq "UNC_Homo" or $genotype eq "P_Homo_var"){
		return 1;
	} else {
		return 0;
	}
}

sub exclude(){
	my ($line,$config) = @_;
	my @sample_ids = $config->get_sample_ids();
	foreach my $sample (@samples){
		my @cols = split(/[$sep]/,$line);
		if (&IsVariant($cols[$config->get_sample_parameter($sample,"Genotype")])){
			return 1;
		}
	}
	return 0;
}

exit;