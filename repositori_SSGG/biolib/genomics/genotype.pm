#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with bioinformatics analysis information
# @Author: Arbol
# @Contributors: JM Rosa
########################################################

package genotype;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Exporter; # Use this library so that scripts invoking functions of this library don't need to include the library name
use Scalar::Util;

# -----------------------
# -- Global variables  --
# -----------------------
our @ISA = qw(Exporter);
our @EXPORT = qw (get_genotype);
our $SNP_af_threshold = 0.12;
our $INDEL_af_threshold = 0.10;

sub get_af_threshold ($) {
	my ($class,$variant) = @_;
	if ($variant eq "SNP"){
		return $SNP_af_threshold;
	} elsif ($variant eq "INDEL"){
		return $INDEL_af_threshold;
	} else {
		die "ERROR: Illegal variant type $variant. Should be either SNP or INDEL!!!\n";
	}
}

sub get_depth_threshold () {
	my ($class) = @_;
	return 3;
}

sub get_genotype ($$$$) {
	my ($class,$depth,$freq,$ref,$var) = @_;
	
	# Normalize allele frequency
	my $ratio;
	if (Scalar::Util::looks_like_number($freq)){
		$ratio = sprintf('%.2f', $freq);
	} else {
		$ratio = 0.00;
	}
	
	# Normalize depth value
	my $rev_depth;
	if (Scalar::Util::looks_like_number($depth)){
		$rev_depth = $depth;
	} else {
		$rev_depth = 0;
	}
	
	my $value;
	
	# If there is no coverage for the position, simply report an "Uncovered" phenotype
	if ($rev_depth == 0){
		$value = "Uncovered,-"
	} else {
		# Different criteria have been stablished for indels and SNVs:
		if (length($ref) == 1 and length($var) == 1){
			# SNVs:
			if ($ratio == 0) {
				$value = "Homo_ref,0";
			}
			elsif ($ratio > 0 && $ratio <= 0.12) {
				$value = "P_Homo_ref,$ratio";
			}
			elsif ($ratio > 0.12 && $ratio < 0.35) {
				$value = "UNC_Hetero,$ratio";
			}
			elsif ($ratio >= 0.35 && $ratio < 0.65) {
				$value = "P_Hetero,$ratio";
			}
			elsif ($ratio >= 0.65 && $ratio < 0.85) {
				$value = "UNC_Homo,$ratio";
			}
			else {
				$value = "P_Homo_var,$ratio";
			}
		} else {
			# Indels:
			if ($ratio == 0) {
				$value = "Homo_ref,0";
			}
			elsif ($ratio > 0 && $ratio <= 0.3) {
				$value = "UNC_Hetero,$ratio";
			}
			elsif ($ratio > 0.3 && $ratio < 0.6) {
				$value = "P_Hetero,$ratio";
			}
			else {
				$value = "P_Homo_var,$ratio";
			}
		}
	}
	
	return $value;
}

sub IsVariant(){
	my ($class,$genotype) = @_;
	if ($genotype eq "UNC_Hetero" or $genotype eq "P_Hetero" or $genotype eq "UNC_Homo" or $genotype eq "P_Homo_var"){
		return 1;
	} else {
		return 0;
	}
}

sub Zygosity()
# Function returning 0 for homozygous for the reference, 1 for heterozygous and 2 for homozygous for the variant
{
	my ($class,$genotype) = @_;
	if ($genotype eq "UNC_Hetero" or $genotype eq "P_Hetero"){
		return 1;
	} elsif ($genotype eq "UNC_Homo" or $genotype eq "P_Homo_var"){
		return 2;
	} else {
		return 0;
	}
}

sub Chromosomes()
# Function that returns the human chromosomes names
{
	my @chroms = (1..22, "X", "Y", "M");
	return \@chroms;
}

sub GetChromosomeNumber($$)
# Function that returns a sorted number for each chromosome (X=23, Y=24, M=25)
{
	my ($class,$chrom) = @_;
	if ($chrom =~ m/chr/i){
		$chrom = substr($chrom,3,length($chrom)-2);
	}
	if ($chrom eq "X" or $chrom eq "x"){
		return 23;
	} elsif ($chrom eq "Y" or $chrom eq "y"){
		return 24;
	} elsif ($chrom eq "M" or $chrom eq "m"){
		return 25;
	}
	else {
		return int($chrom);
	}
}

1;