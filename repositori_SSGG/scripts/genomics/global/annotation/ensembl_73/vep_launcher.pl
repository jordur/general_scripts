#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Variant Effect Predictor launcher for VEP 73,
# @Version: Ensembl 73 launcher v1.0
# @Status: Develop environment
# @Author: gmarco
# @Last modification: 15/10/13
# @Contributors: Arbol 
########################################################

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(firstidx);
use Cwd;
use File::Which;

my $var_predictor = 'variant_effect_predictor.pl';
my $path_to_annot = dirname(which($var_predictor));
my $vep_config = $path_to_annot.'/vep_config/vep_template.config';
my $vep_plugins = $path_to_annot.'/vep_config/Plugins';
my $condel_plugin = "Condel,".$vep_plugins."/config/Condel/config,b";
#my $vep_config = $ENV{'VEP_CONFIG'};

main();

sub main {
	my %options=();
	my ($vcf_input,$vcf_gz,$vcf_tabix,$fields, $samples, $output, $threads);
	my $buffer_size = 25000;
	
	getopts("i:o:t:b:", \%options);
	
	print "\n=================================================================
	Variant Effect Predictor Launcher\n";
	if (not defined $options{i} or not defined $options{o} or not defined $options{t}){
		print "\tYou must specify all tags!n\n";
		print "-i *.vcf file ABSOLUTE PATH ! ie: /input/variants.vcf\n";
		print "-o output_file ABSOLUTE PATH ! ie: /results/output.vcf\n";
		print "-t number_of_threads ie: -t 4\n";
		print "-b buffer size ie: -b 30000. Default value: 25000 variants loaded in memory.\n";
		print "=================================================================\n\n";
		exit;
	}
	
	if (defined $options{b}){		
		#Check that the buffer size is a numeric value
		die "ERROR: The buffer size must be a numeric value !!!\n\n" if not looks_like_number($options{b});
		$buffer_size = $options{b};		
	} 

	print "=================================================================\n";

	#Asign vcf adn gz files to vars to make script launch easier to read.
	$vcf_input = $options{i};
	$vcf_gz = $vcf_input.".gz";
	$vcf_tabix = $vcf_gz.".tbi";
	$output = $options{o};
	$threads = $options{t};

	#Check input file exists
	if (not -e $options{i}){
		die "ERROR: Couldn't read VCF input file.\n";
	}
	
	#Check if output directory exists
	my $output_dir = dirname($options{o});
	if (not -d $output_dir){
		print "WARNING: Couldn't find output directory ! \nCheck that the output path: $output_dir is correct !!\n";
		print "Creating output directory...\n";
		system ("mkdir -p $output_dir");
	}
	
	#Remove $options{i}.gz file.
	unlink("$options{i}.gz");
	
	#Remove Tabix index file
	unlink("$options{i}.gz.tbi");
	
	#Create bgzip file
	print `bgzip -c $options{i} > $options{i}.gz`;
	
	#Index file with Tabix
	print `tabix -p vcf $options{i}.gz`;
	
	#Check gz file
	die "ERROR: Couldn't find gz file !" if !-e $vcf_gz;

	#Check tabix tbi file	
	die "ERROR: Couldn't find tabix tbi file ! Perhaps your VCF input file is not sorted." if !-e $vcf_tabix;
	
	#Check an output path is given
	die "ERROR: You must specify the full output path including file name ie: -o /home/user/output/annotation.vep\n\n" if not defined $options{o};
	
	#Check that a number of threads is specified
	die "ERROR: You must specify the number of threads. ie: -t 4\n\n" if not defined $options{t};

	#Check that the number threads is a numeric value
	die "ERROR: The number of threads must be a numeric value !!!\n\n" if not looks_like_number($options{t});
	
	#We call the VEP launcher with vcf_input, vcf_gz, output
	launch_vep($vcf_input,$output,$vcf_gz,$vcf_tabix,$samples,$threads, $buffer_size);

}

sub launch_vep {
	my ($vcf_input,$output,$vcf_gz,$vcf_tabix,$samples,$threads, $buffer_size) = @_;
	my $output_summary = $output."_summary.html";
	
	#We generate the output column order of the VEP result file.
	my $fields = "SYMBOL,Gene,Gene_description,HGVSc,ALLELE_NUM,INTRON,EXON,HGVSp,Consequence,GMAF,AFR_MAF,AMR_MAF,ASN_MAF,EUR_MAF,InterPro_IDs,InterPro_descriptions,HGMD_info,Related_publication,Existing_variation,Feature_type,Feature,RefSeq,CCDS,CANONICAL,Conservation,Grantham_distance,Condel,SIFT,PolyPhen,DOMAINS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,Flanking_sequence";
	
	print "Launching VEP 73 annotation:\n";
	print "VEP start on: ",`date`,"\n";
	print "Input: $vcf_input\n";
	print "Output: $output\n";
	print "Threads: $threads\n";
	print "Buffer size: $buffer_size\n";
	print "=================================================================\n";
	
	#We launch VEP
	system("$var_predictor -i $vcf_input -o $output -dir_plugins $vep_plugins -plugin $condel_plugin -config $vep_config -plugin Biobase,$vcf_gz -fields $fields -fork $threads -buffer_size $buffer_size");
	print "\n$var_predictor -i $vcf_input -o $output -dir_plugins $vep_plugins -plugin $condel_plugin -config $vep_config -plugin Biobase,$vcf_gz -fields $fields -fork $threads -buffer_size $buffer_size";
	
	print "=================================================================\n";	
	print "DONE: VEP End time: ",`date`,"\n";	
	
	#Remove trash
	#`rm $output_summary`;
	`rm $vcf_gz $vcf_tabix`;
		
	print "Removed temporary files.\n";
	print "DONE: VEP parsing output.\n";
	print "Final annotated file: $output\n";
	print "=================================================================\n"; 
}
