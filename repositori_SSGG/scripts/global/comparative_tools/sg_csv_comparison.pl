#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Compare CSV from different samples
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
sub analysing_field_8 ($$);
sub get_SB ($$$);

my $input_id;
my $output_id;

GetOptions (	
				"l=s"	=> \$input_id,	
				"o=s"	=> \$output_id,	
);

print "\n\n=================================================================
Getting Vars from different samples in csv format script v 1.0\n";

if (not defined $input_id) {
	die "
			Options:
		-l Configuration file (name csv_file_path)
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;
 
if (not defined $output_id) {
	$output_id = "comparison";
}
 
open (LOG, ">$output_id.log") or die "Can´t create log file\n";
print LOG "sg_csv_comparison.pl -l $input_id -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n"; 

open (LIST, "< $input_id") or die "Can´t open $input_id\n";
my %columns;

my %Variants;

while (my $line = <LIST>) {
	chomp $line;
	my ($name, $file) = split (/\s+/, $line);
	open (VARS, "< $file") or die "Can't open $file\n";
	while (my $var = <VARS>) {
		if ($var =~ m/^Gene_name/ or $var =~ m/;;;/) { next;}
		my (@data1, $desc);
		if ($var =~ m/\"/) {
			my @data = split (/\"/, $var);
			@data1 = split (/\;/, $data[0]);
			push (@data1, $data[1]);
			$desc = '"'.$data[1].'"';
			if (not defined $data[2]) { print "$var\n"; }
			my @data2 = split (/\;/, $data[2]);
			push (@data1, @data2);
		}
		else {
			@data1 = split (/\;/, $var);
			$desc = '-';
		}
		my $transcript = 'TRANSCRIPT';
		if (defined $data1[11]) { $transcript = $data1[11]; }
		
		if (defined $Variants{$data1[1]}{$data1[2]}{$data1[3]}{$data1[4]}{$transcript}) {
			$Variants{$data1[1]}{$data1[2]}{$data1[3]}{$data1[4]}{$transcript}{samples} .= ",$name";
		}
		else {
			$Variants{$data1[1]}{$data1[2]}{$data1[3]}{$data1[4]}{$transcript}{samples} = $name;
			my $string = join (";", $data1[0],$data1[1],$data1[2],$data1[3],$data1[4],$data1[8],$data1[9],$data1[10],$data1[11],$desc,$data1[14],$data1[15],$data1[16], $data1[17],$data1[17],$data1[19],$data1[20],$data1[24],$data1[25]);
			$Variants{$data1[1]}{$data1[2]}{$data1[3]}{$data1[4]}{$transcript}{string} = $string;
		}
	}
	close VARS;
} 
close LIST;

open (OUT, ">$output_id.csv") or die "Can´t create $output_id.csv\n";
print OUT "#Gene_name;Chr;Position;Reference_genotype;Sample_genotype;Variant_ID;Biobase_ID;Ensembl_Gene_ID;Ensembl_Transcript_ID;Gene_description;Consequence_on_transcripts;CDS_position;aa_position;aa_change;GERP_conservarion_score;HGVs_ID;Condel_Prediction;Interpro_ID;Interpro_Description;Samples\n";

foreach my $chr (keys %Variants) {
	foreach my $pos (keys %{$Variants{$chr}}) {
		foreach my $ref (keys %{$Variants{$chr}{$pos}}) {
			foreach my $var (keys %{$Variants{$chr}{$pos}{$ref}}) {
				foreach my $transcript (keys %{$Variants{$chr}{$pos}{$ref}{$var}}) {
					print OUT "$Variants{$chr}{$pos}{$ref}{$var}{$transcript}{string};$Variants{$chr}{$pos}{$ref}{$var}{$transcript}{samples}\n";
				}
			}
		}
	}
}

close OUT;

print LOG "\nProcess finished... ".localtime()."\n";
close LOG;

print  "\nProcess finished... ".localtime()."\n";
exit;
