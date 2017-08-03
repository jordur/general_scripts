#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to Filter Variants
# @Author: JM Rosa 
# @Contributors: Arbol
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
use FileHandle;
use analysis_info;
use Scalar::Util;

my $list_id;
my $output_id;

GetOptions (	
				"l=s"	=> \$list_id,	
				"o=s"	=> \$output_id,	
);

print "\n\n=================================================================
Script for comparing variants from tumour samples · v 1.1\n";

if (not defined $list_id) {
	die "
			Options:
		-l File containing config info
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

#if (not defined $output_id){get_output_name($list_id);}

&main ();

sub main () {
	my %sample_info;
	#Getting sample information from config file
	my $n_types;
	my $n_samples = &get_samples($list_id, \%sample_info, \$n_types);
	my @samples = split (/,/, $n_samples);
	my @types = split (/,/, $n_types);
	my %genes;
	my %variants;
	foreach my $type (@types){
		foreach my $sample (@samples){
			print "Analising variants from $sample_info{$sample}{files}{$type}{name}\n";
			get_data ($type, $sample, \%sample_info, \%variants, \%genes);
		}
		print "Printing results ${output_id}_$type\n";
		print_output ($type, \@samples, \%sample_info, \%variants, \%genes);
	}
}

print  "\nProcess finished... ".localtime()."\n";

sub print_output {
	my ($type, $samples, $sample_info, $variants, $genes) = @_;
	my $gene_output = &get_out_file_handle("${output_id}_${type}_gene_table");
	print $gene_output &print_gene_header($samples);
	my $var_output = &get_out_file_handle("${output_id}_${type}_variant_table");
	print $var_output &print_var_header($samples);
	#Looking into the genes
	foreach my $gene (sort {$a cmp $b} keys %{$genes}){
		if (defined $$genes{$gene}{$type}{value} and $$genes{$gene}{$type}{value} == 1){
			#Printing Gene summary
			print $gene_output "$gene\t$$genes{$gene}{description}";
			foreach my $sample (@{$samples}){
				my $number = ($$genes{$gene}{$type}{samples}{$sample}{variants} || "0");
				print $gene_output "\t$number";
			}
			print $gene_output "\n";
			#Printing Variant Summary
			foreach my $pos (sort {$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$a} <=> $$variants{$gene}{$type}{$$genes{$gene}{chr}}{$b}} keys %{$$variants{$gene}{$type}{$$genes{$gene}{chr}}}) {
				foreach my $var (sort {$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$a} <=> $$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$b}} keys %{$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}}) {
					$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$var}{consequences} =~ s/,$//;
					print $var_output "$gene\t$$genes{$gene}{description}\t$$genes{$gene}{chr}\t$pos\t$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$var}{ref}\t$var\t$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$var}{existing}\t$$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$var}{consequences}";
					foreach my $sample (@{$samples}){
						my $type = ($$variants{$gene}{$type}{$$genes{$gene}{chr}}{$pos}{$var}{samples}{$sample}{type} || "-");
						print $var_output "\t$type";
					}
				print $var_output "\n";
				}
			}
		}
		else{
			next;
		}
	}
	return 1;
}


sub get_data {
	my ($type, $s_name, $sample_info, $variants, $genes) = @_;
	my $config = analysis_info_tumours->get_config(psv_file => $$sample_info{$s_name}{files}{$type}{name});
	my $input = &get_in_file_handle($$sample_info{$s_name}{files}{$type}{name});
	while (my $line = <$input>){
		chomp $line;
		if ($line =~ m/^#/){next;}
		my $sep = $$sample_info{$s_name}{files}{$type}{sep};
		my @data = split (/[$sep]/, $line);
		my $gene = $data[$config->{HGNC_symbol}];
		my $chr = $data[$config->{chr}];
		my $position = $data[$config->{position}];
		my $ref = $data[$config->{ref_allele}];
		my $var = $data[$config->{var_allele}];
		my $existing = $data[$config->{existing_col}];
		my $var_type = $data[$config->{var_type}];
		my $consequence = $data[$config->{consequence_col}];
		my $description = $data[$config->{gene_description}];
		
		if (defined $$variants{$gene}{$type}{$chr}{$position}{$var}{samples}{$s_name}{type}){
			if (not defined $$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} or $$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} !~ m/$consequence/){ 
				$$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} .= "$consequence,";
			}
		}
		else {
			$$variants{$gene}{$type}{$chr}{$position}{$var}{samples}{$s_name}{type} = $var_type;
			$$variants{$gene}{$type}{$chr}{$position}{$var}{ref} = $ref;
			$$variants{$gene}{$type}{$chr}{$position}{$var}{existing} = $existing;
			
			if (not defined $$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} or $$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} !~ m/$consequence/){ 
				$$variants{$gene}{$type}{$chr}{$position}{$var}{consequences} .= "$consequence,";
			}
			
			if (defined $$genes{$gene}{$type}){
				if (defined $$genes{$gene}{$type}{samples}{$s_name}{variants}){
					$$genes{$gene}{$type}{samples}{$s_name}{variants} += 1;
				}
				else {
					$$genes{$gene}{$type}{samples}{$s_name}{variants} = 1;
				}
			}
			else {
				$$genes{$gene}{$type}{samples}{$s_name}{variants} = 1;
				$$genes{$gene}{description} = $description;
				$$genes{$gene}{chr} = $chr;
				$$genes{$gene}{$type}{value} = 1;
			}
		}
	}
	return 1;
}

sub get_samples {
	my ($list, $samples, $types) = @_;
	my $input = &get_in_file_handle($list);
	my %tmp_samples;
	my %tmp_types;
	while (my $line = <$input>) {
		chomp $line;
		my @data = split (/\s+/, $line);
		if (defined $$samples{$data[0]}{files}{$data[1]}) {die "File $data[2] already exists: $$samples{files}{$data[1]}\n";}
		else {
			$tmp_samples{$data[0]} = 1;
			$$samples{$data[0]}{files}{$data[1]}{name} = $data[2];
			$$samples{$data[0]}{output}{$data[1]}{output} = $data[0].'_'.$data[1];
			$$samples{$data[0]}{files}{$data[1]}{sep} = $data[3];
		}
		if (not defined $tmp_types{$data[1]}) {
			$tmp_types{$data[1]} = 1;
		}
	}

	foreach my $type (sort {$a cmp $b} keys %tmp_types) {
		$$types .= "$type,";
	}
	$$types =~ s/,$//;
	
	my $results;
	foreach my $sample (sort {$a cmp $b} keys %tmp_samples) {
		$results .= "$sample,";
	}
	$results =~ s/,$//;
	return $results;
}

sub print_gene_header () {
	my $samples = shift;
	my $line = "#Gene_id\tGene_desc";
	foreach my $sample (@{$samples}) {
		$line .= "\t$sample"
	}
	$line .= "\n";
	return $line;
}

sub print_var_header () {
	my $samples = shift;
	my $line = "#Gene_id\tGene_desc\tChr\tPosition\tRef\tVar\tExisting_var\tConsequence";
	foreach my $sample (@{$samples}) {
		$line .= "\t$sample"
	}
	$line .= "\n";
	return $line;
}

sub get_out_file_handle {
	my $file = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open("> $file.csv") or die("ERROR: Could not write to output file $file.csv\n");

	return $out_file_handle;
}

sub get_in_file_handle {
	
	my $file = shift;
	
    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($file)) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $file, "\n") unless -e $file;
        
        $in_file_handle->open( $file ) or die("ERROR: Could not read from input file ", $file, "\n");
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...");
    }
    
    return $in_file_handle;
}

exit;
