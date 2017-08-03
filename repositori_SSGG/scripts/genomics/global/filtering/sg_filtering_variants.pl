#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Filter Variants
# @Author: Arbol & Rosa
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
use variants_info;
use file_handle;
use Scalar::Util;
use Number::FormatEng qw(:all);
use genomics::genotype;

my $psv_id;
my $sep_id = '|';
my $output_id = "filtered";
my $type = "both";
my $depth = 2;
my $freq_snvs = genotype->get_af_threshold("SNP");
my $freq_indels = genotype->get_af_threshold("INDEL");
my $sb = 0;
my $bioscope_sb_up = 1e+100;
my $bioscope_sb_low = 0;
my $conseq_unknown = 0;
my $conseq_known = 0;
my $only_hgmd = 0;
my $not_annot = 0;
my $all_maf = 0.01;

GetOptions (
				"c=s"	=> \$psv_id,	
				"d=s"	=> \$sep_id,
				"o=s"	=> \$output_id,	
				"t=s"	=> \$type,	
				"depth=f"	=> \$depth,
				"freq_snvs=f"	=> \$freq_snvs,
				"freq_indels=f"	=> \$freq_indels,
				"sb=f"	=> \$sb,
				"bioscope_up=f"	=> \$bioscope_sb_up,
				"bioscope_low=f"	=> \$bioscope_sb_low,
				"conseq_unknown"	=> \$conseq_unknown,
				"conseq_known"	=> \$conseq_known,
				"only_hgmd"	=> \$only_hgmd,
				"not_annot"		=> \$not_annot,
				"maf=f"	=> \$all_maf,
);

print "\n\n=================================================================
Script for filtering out variants from SG output variants files (both annotated -psv- and not annotated -vcf- files supported) :: v 1.4\n";

if (not defined $psv_id) {
	die "
		Options:
	-c file_name        Input file name (both SG annotated files -psv- or variants files -vcf- are supported)
	    For not annotated files:
	-not_annot          Option for filtering not annotated variants files (vcf files, default psv)
	-bioscope_up value  Maximum BioScope strand bias value (by default, only values over 1e+100 will be filtered out)
	-bioscope_low value Minimum BioScope strand bias value (0 by default, that is, nothing will be filtered out)
		For annotated files:
	-d separator        Field separator (default = '|')
	-t type             Type of variant (known, unknown or both, default = both)
	-o output_name      Output name (optional, default name is filtered_)
	-conseq_known       Option for filtering known variants after consequence types (only ESSENTIAL_SPLICE_SITE, FRAMESHIFT_CODING, SPLICE_SITE, CODING_UNKNOWN, NON_SYNONYMOUS_CODING, PARTIAL_CODON, STOP_GAINED, STOP_LOST, CODING_SEQUENCE_VARIANT, FEATURE_ELONGATION, FEATURE_TRUNCATION, SPLICE_REGION_VARIANT, STOP_RETAINED_VARIANT remain, removing all WITHIN_NON_CODING_GENE and NMD_TRANSCRIPT). Only in case of annotated variants available (default is 0=FALSE)
	-conseq_unknown     Option for filtering unknown variants after consequence types (only ESSENTIAL_SPLICE_SITE, FRAMESHIFT_CODING, SPLICE_SITE, CODING_UNKNOWN, NON_SYNONYMOUS_CODING, PARTIAL_CODON, STOP_GAINED, STOP_LOST, CODING_SEQUENCE_VARIANT, FEATURE_ELONGATION, FEATURE_TRUNCATION, SPLICE_REGION_VARIANT, STOP_RETAINED_VARIANT remain, removing all WITHIN_NON_CODING_GENE and NMD_TRANSCRIPT). Only in case of annotated variants available (default is 0=FALSE)
	-only_hgmd          Option for only including, apart from unknown variants, those HGMD described known variants. Only in case of annotated variants available
	-maf per_unit_value Filters out all variants whose defined ALL_MAF population frequencies lie above the \"maf\" per unit expressed value (default >=0.01)
		For both annotated or not annotated files:
	-depth value        Minimum coverage value (default > 2)
	-freq_snvs value    Minimum SNPs allele freq value (default > 0.12)
	-freq_indels value  Minimum indels allele freq value (default > 0.10)
	-sb value           Minimum Strand bias value in Phred scale (-10*log10(p)) for GATK and Varscan (by default, a value of 0 is taken, so that no variant will be filtered out for this reason)
\n".localtime()."\n=================================================================\n\n";
}

&main ();
print  "\nProcess finished... ".localtime()."\n";

sub main () {
	# For allowing input parameters in scientific format...
	$sb = format_eng($sb);
	$bioscope_sb_up = format_eng($bioscope_sb_up);
	$bioscope_sb_low = format_eng($bioscope_sb_low);
	$all_maf = format_eng($all_maf);

	# Config will be loaded, ouput file handles are created
	my $config;
	my %output;
	if ($not_annot){
		$config = variants_info->get_config(vcf_file => $psv_id);
		%output = &get_output($config,"vcf");
	}
	else {
		$config = analysis_info->get_config(type => $type, psv_file => $psv_id);
		%output = &get_output($config, "psv");
	}
	my $input = &get_in_file_handle($psv_id);
	my $string;
	my $tip;
	while (my $line = <$input>){
		chomp $line;
		if ($line =~ m/^#/) {
			foreach(keys(%output)){
				my $output = $output{$_};
				print $output $line."\n";
			}
		} else {
			$string = get_results($line, $config, \$tip);
			my $output = $output{$tip};
			print $output "$string";		
		}
	}
}

sub get_sufix() {
	my ($number,$default) = @_;
	if ($number == 1){ return 'strand_bias';}
	elsif ($number == 2){ return 'low_depth';}
	elsif ($number == 3){ return 'low_freq';}
	elsif ($number == 4){ return 'consequence';}
	elsif ($number == 5){ return 'non_hgmd_vars';}
	elsif ($number == 6){ return 'MAF_polymorphisms';}
	else {return $default;}
}

sub conseq_filter($){
	my ($data) = @_;
	# Coyote's filtering criteria:
	#if ($data->{Variant_effect} =~ m/NMD_TRANS/g) {return 4;}
	#if ($data->{Variant_effect} =~ m/WITHIN_NON/g) {return 4;}
	#if ($data->{Variant_effect} eq 'INTRONIC') {return 4;}
	#if ($data->{Variant_effect} eq 'DOWNSTREAM') {return 4;}
	#if ($data->{Variant_effect} eq 'UPSTREAM') {return 4;}
	# New filtering criteria:
	my $variant_effect = $data->get_parameter("Variant_effect");
	if ($variant_effect =~ m/WITHIN_NON_CODING_GENE/g) {return 1;}
	if ($variant_effect =~ m/NMD_TRANSCRIPT/g) {return 1;}
	if (not($variant_effect =~ m/splice_donor_variant/g) 
		and not($variant_effect =~ m/splice_acceptor_variant/g) 
		and not($variant_effect =~ m/stop_gained/g) 
		and not($variant_effect =~ m/frameshift_variant/g) 
		and not($variant_effect =~ m/stop_lost/g) 
		and not($variant_effect =~ m/inframe_insertion/g) 
		and not($variant_effect =~ m/inframe_deletion/g) 
		and not($variant_effect =~ m/splice_region_variant/g) 
		and not($variant_effect =~ m/incomplete_terminal_codon_variant/g) 
		and not($variant_effect =~ m/missense_variant/g) 
		and not($variant_effect =~ m/initiator_codon_variant/g)) {return 1;}
	return 0;
}

sub get_results ($$$) {
	my ($line, $config, $tip) = @_;
	
	# Get results from line of either vcf or psv file
	my $results;
	if ($not_annot){
		$results = variants_info->get_fields_from_vcf($line,$config);
		$$tip = 'not_annotated';
		my $number;
		$number = &filtering ($config, $results);
		$$tip = &get_sufix($number,$$tip);
		return $results->put_fields_into_vcf($config);
	} else {
		$results = analysis_info->get_fields_from_psv($line,$config);
		$$tip = &get_type($config, $results);
		my $number; 
		
		if ($config->{type} eq 'known') {
			if ($$tip eq 'unknown') {
				$$tip = 'failed';
			}
			else {
				$number = &filtering ($config, $results);
				$$tip = &get_sufix($number,$$tip);
			}
		}
		elsif ($config-> {type} eq 'unknown') {
			if ($$tip eq 'known') {
				$$tip = 'failed';		}
			else {
				$number = &filtering ($config, $results);
				$$tip = &get_sufix($number,$$tip);
			}
		}
		else {
			$number = &filtering ($config, $results);
			$$tip = &get_sufix($number,$$tip);
		}
		return $results->put_fields_into_psv($config);
	}
}

sub filtering () {
	my ($config, $data)	= @_;
	my $result = 0; # 0=not being filtered, 1=filtered by strand bias, 2=filtered by depth, 3=filtered by freq
	# 4=filtered by consequence, 5=filtered for not being an hgmd described known mutation,
	# 6=filtered by population freq
	
	# Normalize strand-bias values:
	my $SB_normalized;
	my $SB = $data->get_parameter("Strand_Bias");
	if ((Scalar::Util::looks_like_number($SB)?"true":"false") eq 'false'){
		if ($SB eq "-"){
			$SB_normalized = 0;
		} elsif ($SB eq "+++"){
			$SB_normalized = 1e+99;
		} elsif ($SB eq "---"){
			$SB_normalized = 1e-99;
		} else {
			print "WARNING: $SB is NOT numeric!!\n";
		}
	}
	else {
		$SB_normalized = $SB;
	}
	
	if ($not_annot){
		# In not annotated files, different strand bias criteria can be applied for BioScope variants
		# GATK strand bias filtering
		if ($data->get_parameter("gatk") ne "0"){
			if ($SB_normalized < $sb){
				return 1;
			}
		}
		# BioScope strand bias filtering
		else {
			if ($SB_normalized > $bioscope_sb_up or $SB_normalized < $bioscope_sb_low){
				return 1;
			}
		}
	} else {
		if ($SB_normalized < $sb){
			return 1;
		}
	}
	
	# Samples related filters (now only depth and freq)
	my $sum_depth = 0;
	my $sum_freq = 0;
	my @sample_ids = $config->get_sample_ids();
	foreach my $sample_name (@sample_ids){
		
		# Depth filtering
		if ($data->get_sample_parameter($sample_name, "Depth") ne '-' and $data->get_sample_parameter($sample_name, "Depth") > $depth) {$sum_depth++;}

		# Freq filtering
		# Different filtering will be carried out depending on wether the variant is snp or indel
		if (length($data->get_parameter("Ref_Allele")) == 1 and length($data->get_parameter("Var_Allele")) == 1){
			if ($data->get_sample_parameter($sample_name, "Freq") ne '-' and $data->get_sample_parameter($sample_name, "Freq") > $freq_snvs) {
				$sum_freq++;
			}
		} else {
			if ($data->get_sample_parameter($sample_name, "Freq") ne '-' and $data->get_sample_parameter($sample_name, "Freq") > $freq_indels) {
				$sum_freq++;
			}
		}
	}
	if ($sum_depth == 0) {
		return 2;
	}
	if ($sum_freq == 0) {
		return 3;
	}
	
	# Annotation-dependent filtering
	if (not $not_annot){
		
		# Known and unknown variants filtering
		if ($data->get_parameter("Existing_variation") eq "-" and $data->get_parameter("HGMD_info") eq "-"){
			# Filtering after consequence only in case of UNKNOWN VARIANTS
			if ($conseq_unknown){
				if (conseq_filter($data)) {return 4;}
			}
		} else {
			if ($conseq_known){
				if (conseq_filter($data)) {return 4;}
			}
			if ($only_hgmd){
				# Filtering of KNOWN VARIANTS (only HGMD described variants will remain)
				if ($data->{HGMD_info} eq "-") {return 5;}
			}
		}
		
		# Populations frequencies filtering
		if (Scalar::Util::looks_like_number($data->get_parameter("ALL_MAF"))){
			if ($data->get_parameter("ALL_MAF") >= $all_maf) {return 6;}
		}
	}
	return $result;
}

sub get_type () {
	my ($config, $data) = @_;
	if ($type ne "both"){
		if ($data->{existing_variation} ne "-"){ #Former value: =~ m/^rs/){
			return 'known';
		}
		else {
			return 'unknown';
		}
	}
	return $type;
}

sub get_output () {
	my ( $config, $file_type) = @_;

	my %output;
	
	if ($not_annot){
		$output{not_annotated} = &get_out_file_handle("${output_id}_not_annotated", $file_type);
	} else {
		if ($type ne "both"){
			if ($config->{print_known} eq 'true'){
				$output{known} = &get_out_file_handle("${output_id}_known", $file_type);
			} 
			if ($config->{print_unknown} eq 'true'){
				$output{unknown} = &get_out_file_handle("${output_id}_unknown", $file_type);
			}
		} else {
			$output{both} = &get_out_file_handle("${output_id}_both", $file_type);
		}
		$output{consequence} = &get_out_file_handle("${output_id}_consequence", $file_type);
		$output{non_hgmd_vars} = &get_out_file_handle("${output_id}_non_hgmd_vars", $file_type);
		$output{MAF_polymorphisms} = &get_out_file_handle("${output_id}_MAF_polymorphisms", $file_type);
	} 
	
	$output{low_depth} = &get_out_file_handle("${output_id}_low_depth", $file_type);
	$output{strand_bias} = &get_out_file_handle("${output_id}_strand_bias", $file_type);
	$output{low_freq} = &get_out_file_handle("${output_id}_low_freq", $file_type);
	$output{failed} = &get_out_file_handle("${output_id}_failed", $file_type);
	
	return %output;
}

exit;
