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

my $vep_path = $ENV{'VEP73_SCRIPT'};
my $vep_config = $ENV{'VEP73_CONFIG'};

main();

sub main {
	my %options=();
	my ($vcf_input,$vcf_gz,$vcf_tabix,$fields, $samples, $output, $threads);
	my $buffer_size = 25000;
	
	getopts("i:o:t:b:", \%options);
	
	print "\n=================================================================
	Variant Effect Predictor Launcher 73\n";
	if (not defined $options{i} or not defined $options{o} or not defined $options{t}){
		print "\tYou must specify all tags!n\n";
		print "-i *.vcf file ABSOLUTE PATH ! ie: /input/variants.vcf\n";
		print "-o output_file ABSOLUTE PATH ! ie: /results/output.psv\n";
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
	
	#We obtain the samples from VCF input file
	$samples = &get_samples_col($options{i});
	
	if ($samples eq "") {
		die "Samples header line is missing. Input VCF file format must be correct !\n";
	}
	
	#We call the VEP launcher with vcf_input, vcf_gz, output
	launch_vep($vcf_input,$output,$vcf_gz,$vcf_tabix,$samples,$threads, $buffer_size);

}

sub get_samples_col {
	my $vcf_input = $_[0];
	my $sample_cols="";
	my (@names,@sample_space, @sample_comma, @sample_equal, $strand_bias);
		
	#Open VCF file
	open VCF_INPUT, $vcf_input or die $!;
	
	#Extract the first non header line and get the samples
	#name and the INFO column to count the samples.
	while (<VCF_INPUT>){
		if ($_ =~ /^#Samples/) {
			@sample_space = (split(" ",$_))[1];
		}
	}
	
	close VCF_INPUT;
	
	foreach my $i (@sample_space){
		@sample_comma = split(";",$i);
	}
	
	foreach my $j (@sample_comma){
		push (@sample_equal, (split("=", $j))[1]);
	}

	@names = @sample_equal;
	
	foreach my $sample_name (@names){
		#We generate col names
		$sample_cols .= $sample_name."_Genotype,";
		$sample_cols .= $sample_name."_Depth,";
		$sample_cols .= $sample_name."_Var/Depth,";
	}
	return $sample_cols;
}

sub launch_vep {
	my ($vcf_input,$output,$vcf_gz,$vcf_tabix,$samples,$threads, $buffer_size) = @_;
	my $output_summary = $output."_summary.html";
	
	#We generate the output column order of the VEP result file.
	my $fields = "SYMBOL,Chr,Pos,Ref_Allele,Var_Allele,".$samples."Strand_Bias,Gene,Gene_description,HGVSc,INTRON,EXON,HGVSp,Consequence,GMAF,AFR_MAF,AMR_MAF,ASN_MAF,EUR_MAF,InterPro_IDs,InterPro_descriptions,HGMD_info,Related_publication,Existing_variation,Feature_type,Feature,RefSeq,CCDS,CANONICAL,Conservation,Grantham_distance,Condel,SIFT,PolyPhen,DOMAINS,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,Flanking_sequence";
	
	print "Launching VEP 73 annotation:\n";
	print "VEP start on: ",`date`,"\n";
	print "Input: $vcf_input\n";
	print "Output: $output\n";
	print "Threads: $threads\n";
	print "Buffer size: $buffer_size\n";
	print "=================================================================\n";
	
	#We launch VEP
	system("$vep_path -i $vcf_input -o $output -config $vep_config -plugin vcf_input,$vcf_gz -plugin Biobase,$vcf_gz -fields $fields -fork $threads -buffer_size $buffer_size");
	print "\n$vep_path -i $vcf_input -o $output -config $vep_config -plugin vcf_input,$vcf_gz -plugin Biobase,$vcf_gz -fields $fields -fork $threads -buffer_size $buffer_size\n";
	
	print "=================================================================\n";	
	print "DONE: VEP End time: ",`date`,"\n";	

	#We parse the VEP output
	print "=================================================================\n";
	print "Parsing VEP output...\n"; 	
	vep_output_parser($output, "$output.final");
	
	#Remove trash
	`rm $output $output_summary`;
	`rm $vcf_gz $vcf_tabix`;
	
	#Rename
	`mv $output.final $output`;
	
	print "Removed temporary files.\n";
	print "DONE: VEP parsing output.\n";
	print "Final annotated file: $output\n";
	print "=================================================================\n"; 
}


sub vep_output_parser{

my ($vep_output, $vep_parsed) = @_ ;
my ($gmaf, $gmaf_index, @cols, $cols);

	#Open VEP output file & VEP output parsed file to write.
	open VEP_OUTPUT, "$vep_output";
	open VEP_PARSED, ">$vep_parsed";
	
	#We start parsing VEP output
	while (<VEP_OUTPUT>){
		chomp $_;
		
		#Remove headers
		if ($_ =~ /^(##)/) {next;}

		#Get GMAF column position and print header
		if ($_ =~ /^#/){
			
			@cols = split ("\t", $_);
			$gmaf_index = firstidx { $_ eq 'GMAF' } @cols;
			
			s/SYMBOL/HGNC_symbol/, s/HGVSc/HGVSc_name/, s/HGVSp/HGVSp_name/, s/INTRON/Intron/, 
			s/EXON/Exon/, s/Consequence/Variant_effect/, s/^Feature$/Feature_ID/, s/RefSeq/RefSeq_ID/,
			s/CCDS/CCDS_ID/, s/CANONICAL/Canonical_isoform/, s/Condel/Condel_prediction/, s/GMAF/ALL_MAF/, 
			s/SIFT/SIFT_prediction/, s/PolyPhen/PolyPhen_prediction/, s/DOMAINS/Affected_prot_domains/,
			s/Conservation/Conservation_score/, s/MOTIF_NAME/Regulatory_Motif_name/, 
			s/MOTIF_POS/Regulatory_Motif_position/,s/HIGH_INF_POS/Regulatory_High_Inf_Pos/,
			s/MOTIF_SCORE_CHANGE/Regulatory_Motif_Score_Change/ for @cols;
			
			print VEP_PARSED "",(join "|", @cols), "\n";
			next;
		}
		
		#Parse line
		my @linea = split ("\t",$_);
		chomp @linea;
		
		#If GMAF value is not empty "-" we remove the allele from it.
		if ($linea[$gmaf_index]  ne "-"){ 
			$linea[$gmaf_index] =  (split (":",$linea[$gmaf_index]))[1];
		} 
		print VEP_PARSED "",(join "|", @linea), "\n";
	}
	
	#Close files & handlers
	close VEP_OUTPUT;
	close VEP_PARSED; 
}
