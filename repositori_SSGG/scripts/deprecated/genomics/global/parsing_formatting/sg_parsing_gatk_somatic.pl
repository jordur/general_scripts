#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic GATK in target regions
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
sub get_info ($$$);
sub get_type ($$);
sub get_genotype ($);
sub analysing_field_8 ($$);

my $vcf_id;
my $output_id;
my $normal_id;
my $tumour_id;

GetOptions (	
				"v=s"	=> \$vcf_id,	
				"o=s"	=> \$output_id,	
				"n=s"	=> \$normal_id,	
				"t=s"   => \$tumour_id,			
);

print "\n\n=================================================================
Getting Somatic Vars from VarScan and GATK in target regions script v 1.0\n";

if (not defined $vcf_id or not defined $normal_id or not defined $tumour_id) {
	die "
			Options:
		-v VCF file
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
print LOG "sg_parsing_gatk_somatic.pl -v $vcf_id -n $normal_id -t $tumour_id -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

open (VCF, "< $vcf_id") or die "Can't open $vcf_id\n";
my %columns;
open (OUT, ">$output_id.snvs") or die "Can't create $output_id.snvs\n";
open (ERR, ">$output_id.error") or die "Can't create $output_id.error\n";

print OUT "#Chr\tPosition\tRef_Allele\tVar_Allele\tType\t$tumour_id\_depth\t$tumour_id\_geno\t$tumour_id\_freq\t$normal_id\_depth\t$normal_id\_geno\t$normal_id\_freq\tSB\n";

while (my $line = <VCF>) {
	chomp $line;
	if ($line =~ m/^##/) { next;}
	elsif ($line =~ m/^#CHROM/) {
		$columns{normal} = &getpos ($normal_id, $line);	
		$columns{tumour} = &getpos ($tumour_id, $line);	
	}
	else {
		my @data = split (/\s+/, $line);
		my ($normal_depth, $normal_freq, $tumour_depth, $tumour_freq, $sb, $type, $normal_info, $tumour_info, $normal_geno, $tumour_geno);

		my @vars = split (/\,/, $data[4]);
		my $n = 1;
		foreach my $var (@vars) {
			if ($data[$columns{normal}] ne './.' and $data[$columns{tumour}] ne './.') {
				$normal_depth = &get_info ($data[$columns{normal}], 'depth', $n);
				$normal_info = &get_info ($data[$columns{normal}], 'freq', $n);
				$normal_freq = (split(/\,/, $normal_info))[1];
				$normal_geno = (split(/\,/, $normal_info))[0];
				$tumour_depth = &get_info ($data[$columns{tumour}], 'depth', $n);
				$tumour_info = &get_info ($data[$columns{tumour}], 'freq', $n);
				$tumour_freq = (split(/\,/, $tumour_info))[1];
				$tumour_geno = (split(/\,/, $tumour_info))[0];
				$type = &get_type ($normal_geno, $tumour_geno); 
			}
			elsif ($data[$columns{normal}] ne './.' and $data[$columns{tumour}] eq './.') {
				$normal_depth = &get_info ($data[$columns{normal}], 'depth', $n);
				$normal_info = &get_info ($data[$columns{normal}], 'freq', $n);
				$normal_freq = (split(/\,/, $normal_info))[1];
				$normal_geno = (split(/\,/, $normal_info))[0];
				$tumour_depth = '-';
				$tumour_info = '-,-';
				$tumour_freq = (split(/\,/, $tumour_info))[1];
				$tumour_geno = (split(/\,/, $tumour_info))[0];
				$type = 'unknown'; 
			}
			elsif ($data[$columns{normal}] eq './.' and $data[$columns{tumour}] ne './.') {
				$normal_depth = '-';
				$normal_info = '-,-';
				$normal_freq = (split(/\,/, $normal_info))[1];
				$normal_geno = (split(/\,/, $normal_info))[0];
				$tumour_depth = &get_info ($data[$columns{tumour}], 'depth', $n);
				$tumour_info = &get_info ($data[$columns{tumour}], 'freq', $n);
				$tumour_freq = (split(/\,/, $tumour_info))[1];
				$tumour_geno = (split(/\,/, $tumour_info))[0];
				$type = 'unknown'; 
			}	
			else {
				print ERR "$line\n";
			}
			
			$sb = &analysing_field_8 ($data[7], 'SB=');
			print OUT "$data[0]\t$data[1]\t$data[3]\t$data[4]\t$type\t$tumour_depth\t$tumour_geno\t$tumour_freq\t$normal_depth\t$normal_geno\t$normal_freq\t$sb\n";
		}
		$n++;
	}
} 

close VCF;
close OUT;
close ERR;

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

sub get_info ($$$) {
	my ($line, $match, $n) = @_;
	my @fields =  split (/:/, $line);
	my $ref = (split (/\,/, $fields[1]))[0];
	my $var = (split (/\,/, $fields[1]))[$n];
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
			$value = &get_value ($freq,$ref,$var); 
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
		$value = "UNC_Homo";
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
