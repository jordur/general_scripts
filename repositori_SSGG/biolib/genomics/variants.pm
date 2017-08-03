#!/usr/bin/perl -w
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with genomic variants
# @Author: Arbol
# @Contributors:
########################################################

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;

# -----------------------
# -- Global variables  --
# -----------------------

package variants;
sub new
# Creation of new variant object
{
	my $class = shift; # the first parameter is the class name
	my %args = @_;
	
	my %defaults = (
		chr => 1,
		position => 0,
		ref_allele => '',
		var_allele => '',
		depth => 0,
		genotype => '',
		freq => 0,
		SB => 0,
	);

	my $self = {%defaults,@_};
	
	# Remove the "chr" part from the chromosome name
	$self->{chr} =~ s/^chr//;
	
	# Create and return object
	bless $self, $class;
	return $self;
}

sub create_from_tsv
# Creation of new variant from a tsv line
{
	my ($class, $line) = @_;
	
	my @cols = split ('\t', $line);
	
	# Perform different checks on obtained data to assess its validity
	if ( $#cols < 8 ){
		print "WARNING: Not a valid variant (not enough information for chromosome, position, reference allele, variant allele, depth, genotype, freq, reads+ and reads-)\n";
		return;
	}
	
	my $SB;
	if ($cols[8] == 0){
		$SB = '1e99';
	} else {
		$SB = $cols[7]/$cols[8];
	}
	
	my $self = variants->new(
		chr => $cols[0],
		position => $cols[1],
		ref_allele => $cols[2],
		var_allele => $cols[3],
		depth => $cols[4],
		genotype => $cols[5],
		freq => $cols[6],
		SB => $SB,
    );
    
    # Return object
    return $self;
}

sub compare
# Function that returns "0" if both argument variants are the same, "-1" if variant1 < variant2, and "1" if variant1 > variant2
# Note: two variants are considered the same, if their "chr", "position", "ref_allele" and "var_allele" properties are the same
{
	my ($class, $variant1, $variant2) = @_;
	my $comparison;
	
	my $chrom1 = $variant1->{chr};
	my $chrom2 = $variant2->{chr};
	if ($chrom1 eq "X"){ $chrom1 = 23;}
	elsif ($chrom1 eq "Y"){ $chrom1 = 24;}
	elsif ($chrom1 eq "M"){ $chrom1 = 25;}
	if ($chrom2 eq "X"){ $chrom2 = 23;}
	elsif ($chrom2 eq "Y"){ $chrom2 = 24;}
	elsif ($chrom2 eq "M"){ $chrom2 = 25;}
	
	if ($chrom1 < $chrom2){
		$comparison = -1;
	} elsif ($chrom1 == $chrom2){
		if ($variant1->{position} < $variant2->{position}){
			$comparison = -1;
		} elsif ($variant1->{position} == $variant2->{position}){
			if ($variant1->{ref_allele} lt $variant2->{ref_allele}){
				$comparison = -1;
			} elsif ($variant1->{ref_allele} eq $variant2->{ref_allele}){
				if ($variant1->{var_allele} lt $variant2->{var_allele}){
					$comparison = -1;
				} elsif ($variant1->{ref_allele} eq $variant2->{ref_allele}){
					$comparison = 0;
				} else {
					$comparison = 1;
				}
			} else {
				$comparison = 1;
			}
		} else {
			$comparison = 1;
		}
	} else {
		$comparison = 1;
	}
	
	return $comparison;
}

sub get_first_variant
# Get the first variant of a given array of variants
{
	my ($class, $variants) = @_;
	
	my $first_variant;
	my $i = 0;
	do {
		if (defined($$variants[$i])){
			$first_variant = $$variants[$i];
		} else {
			$i ++;
		}
	} while ((not defined($first_variant)) and $i <= $#{$variants});
	
	# Return "undef" if all elements of array are undef
	if (not defined($first_variant)){
		return undef;
	}
	
	my $pos_first = $i;
	for (my $i=$pos_first;$i<=$#{$variants};$i++){
		if (defined($$variants[$i])){
			if (variants->compare($$variants[$i],$first_variant) == -1){
				$first_variant = $$variants[$i];
				$pos_first = $i;
			}
		}
	}
	    
	# Return object
	return $first_variant;
}

1;