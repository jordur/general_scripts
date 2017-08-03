#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic VarScan in target regions
# @Author: JM Rosa
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

sub getpos ($$);
sub get_info ($$);
sub get_type ($$);
sub get_value ($);
sub analysing_field_8 ($$);
sub get_SB ($$$);

my $indel_id;
my $snp_id;
my $output_id;
my $normal_id;
my $tumour_id;

GetOptions (	
				"i=s"	=> \$indel_id,	
				"s=s"	=> \$snp_id,	
				"o=s"	=> \$output_id,	
				"n=s"	=> \$normal_id,	
				"t=s"   => \$tumour_id,			
);

print "\n\n=================================================================
Getting Somatic Vars from VarScan and GATK in target regions script v 1.0\n";

if (not defined $snp_id or not defined $indel_id or not defined $normal_id or not defined $tumour_id) {
	die "
			Options:
		-s SNPs file
		-i Indels file
		-n Normal sample id
		-t Tumour sample id
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;
 
if (not defined $output_id) {
	$output_id = $normal_id."_".$tumour_id."_gatk";
}
 
open (LOG, ">$output_id.log") or die "Can't create log file\n";
print LOG "sg_parsing_varscan_somatic.pl -s $snp_id -i $indel_id -n $normal_id -t $tumour_id -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

my $command = '$0';
open (INP, "cat $snp_id $indel_id | sed -e 's/^chr//' | sed -e 's/^X/23/' | sed -e 's/^Y/24/' | sed -e 's/^M/25/' | sort -k1n -k2n | sed -e 's/^23/X/' | sed -e 's/^24/Y/' | sed -e 's/^25/M/' | awk '{print \"chr\"$command}' |") or die "Can't open input files\n";
open (OUT, ">$output_id.snvs") or die "Can't create $output_id.snvs\n";
open (ERR, ">$output_id.error") or die "Can't create $output_id.error\n";

print OUT "#Chr\tPosition\tRef_Allele\tVar_Allele\tType\t$tumour_id\_depth\t$tumour_id\_geno\t$tumour_id\_freq\t$tumour_id\_SB\t$normal_id\_depth\t$normal_id\_geno\t$normal_id\_freq\t$normal_id\_SB\n";

# chrom   position        ref     var     normal_reads1   normal_reads2   normal_var_freq normal_gt       tumor_reads1    tumor_reads2    tumor_var_freq  tumor_gt        somatic_status  variant_p_value somatic_p_value tumor_reads1_plus       tumor_reads1_minus      tumor_reads2_plus       tumor_reads2_minus

 while (my $line = <INP>) {
	chomp $line;
	my @data = split (/\s+/, $line);
	my ($normal_depth, $normal_freq, $tumour_depth, $tumour_freq, $tumour_sb, $normal_sb, $type, $normal_info, $tumour_info, $normal_geno, $tumour_geno);
	
	$normal_depth = $data[4] + $data[5];
	$normal_freq = sprintf('%.2f', ($data[5] / $normal_depth));
	$normal_geno = &get_value ($normal_freq);
	$tumour_depth = $data[8] + $data[9];;
	$tumour_freq = sprintf('%.2f', ($data[9] / $tumour_depth));
	$tumour_geno = &get_value ($tumour_freq);
	$normal_sb = '-';
	my $minus_reads = $data[16] + $data[18];
	my $plus_reads = $data[15] + $data[17];
	if ($minus_reads == 0) { 
		$tumour_sb = '+++';
	}
	elsif ($plus_reads == 0) {
		$tumour_sb = '---';
	}
	else {
		$tumour_sb = $plus_reads/$minus_reads;
	}
	
	$type = &get_type ($normal_geno, $tumour_geno); 
	
	my $ref = $data[2];
	my @vars = split(/\//, $data[3]);
	foreach my $var (@vars) {
		if ($var =~ m/-/) { #deletion
			$var = '-';
			$ref = $data[3];
			$ref =~ s/-//;
			$data[1] += 1;
		}
		elsif ($var =~ m/\+/) { #insertion
			$ref = '-';
			$var =~ s/\+//;
			$data[1] += 1;
		}
		print OUT "$data[0]\t$data[1]\t$ref\t$var\t$type\t$tumour_depth\t$tumour_geno\t$tumour_freq\t$tumour_sb\t$normal_depth\t$normal_geno\t$normal_freq\t$normal_sb\n";
	}
 } 

close INP;
close OUT;
close ERR;

sub get_SB ($$$) {
	my ($normal_geno, $tumour_geno, $data) = @_;
	my @info;
	if ($normal_geno =~ m/ref/i) {
		$info[1] = &analysing_field_8 ($data, 'SB=');
		$info[0] = 'N\A';
	}
	elsif ($tumour_geno =~ m/ref/i) {
		$info[0] = &analysing_field_8 ($data, 'SB=');
		$info[1] = 'N\A';
	}
	elsif ($normal_geno eq '-' and $tumour_geno ne '-') {
		$info[1] = &analysing_field_8 ($data, 'SB=');
		$info[0] = 'N\A';
	}
	elsif ($normal_geno ne '-' and $tumour_geno eq '-') {
		$info[0] = &analysing_field_8 ($data, 'SB=');
		$info[1] = 'N\A';
	}
	else {
		$info[0] = &analysing_field_8 ($data, 'SB=');
		$info[1] = &analysing_field_8 ($data, 'SB=');
	}
	return @info
}

sub get_type ($$) {
	my ($normal_gen, $tumour_gen) = @_;
	if ($normal_gen eq $tumour_gen) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/ref/i and $tumour_gen =~ m/ref/i) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/ref/i and $tumour_gen !~ m/ref/i) {
		return "Somatic";
	}
	elsif ($normal_gen !~ m/ref/i and $tumour_gen =~ m/ref/i) {
		return "LOH";
	}
	elsif ($normal_gen !~ m/var/i and $tumour_gen =~ m/var/i) {
		return "LOH";
	}
	elsif ($normal_gen =~ m/hetero/i and $tumour_gen =~ m/hetero/i) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/var/i and $tumour_gen =~ m/hetero/i) {
		return "GOH";
	}
	else {
		return "unknown";
	}
}

sub analysing_field_8 ($$){
	my ($line, $match) = @_;
	my @fields =  split (/;/, $line);
	my $value = '-';
	foreach my $field (@fields) {
		if ($field =~ m/$match/) {
			$field =~ s/$match//;
			$value = $field;
		}
	}
	return $value;
}

sub get_info ($$) {
	my ($line, $match) = @_;
	my @fields =  split (/:/, $line);
	my $ref = (split (/\,/, $fields[1]))[0];
	my $var = (split (/\,/, $fields[1]))[1];
	my $value;

	if ($match eq 'depth') {
		$value = $var + $ref;
	}
	elsif ($match eq 'freq') {
		my $depth = $var + $ref;
		if ($depth == 0) {
			$value = '-,-';
		}
		else {
			my $freq = $var / $depth;
			$value = &get_value ($freq); 
		}
	}
	return $value;
}

sub get_value ($) {
	my $freq = shift;
	my $ratio = sprintf('%.2f', $freq);
	my $value;
	if ($ratio == 0) {
		$value = "Homo_ref";
	}
	elsif ($ratio > 0 && $ratio <= 0.12) {
		$value = "P_Homo_ref";
	}
	elsif ($ratio > 0.12 && $ratio < 0.35) {
		$value = "UNC_Hetero";
	}
	elsif ($ratio >= 0.35 && $ratio < 0.65) {
		$value = "P_Hetero";
	}
	elsif ($ratio >= 0.65 && $ratio < 0.85) {
		$value = "UNC_Hetero";
	}
	else {
		$value = "P_Homo_var";
	}
	return $value;
}

sub getpos ($$) {
	my ($id, $line) = @_;
	my @data = split(/\s+/,$line);
	for (my $i = 0; $i <= scalar(@data); $i++){
		if ($data[$i] eq $id) {
			return $i;
			last;
		}
	}
}

print LOG "\nProcess finished... ".localtime()."\n";
close LOG;

print  "\nProcess finished... ".localtime()."\n";
exit;
