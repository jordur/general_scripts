#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Arbol
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use analysis_info;
use file_handle;
use genomics::genotype;
use File::Path;

# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $input_dominant;
my $input_recessive;
my $sep = '|';
my $sep_out = "\t";
my $outfolder = cwd();
my $prefix = "";
my @zygosity = ("homozygous","heterozygous","compound");
my @inheritance_modes = ("dominant","recessive");
my @file_types = ("genes","variants"); 
my $samples = "";
my $affected = "";
my $zygoresultspath = "zygosity_specific_results";
my $sampleresultspath = "sample_specific_results";
my $compound_dominant = 1; # Number of variants that a sample must at least contain to consider that it's affected in an heterozygous-compound mode for the dominant mode
my $compound_recessive = 2; # Number of variants that a sample must at least contain to consider that it's affected in an heterozygous-compound mode for the recessive mode

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"id=s"	=> \$input_dominant,
				"ir=s"	=> \$input_recessive,
				"fi=s"	=> \$sep,
				"fo=s"  => \$sep_out,
				"s=s"	=> \$samples,
				"a=s"	=> \$affected,
				"o=s"	=> \$outfolder,
				"cd=i"	=> \$compound_dominant,
				"cr=i"	=> \$compound_recessive,
				"p=s"   => \$prefix,
);

print "\n=================================================================
$basename: Script for obtaining candidate genes and variants for dominant and recessive models from SSGG annotation files · v1.1\n";

if (not defined($input_dominant) or $compound_recessive < 2) {
	die "
			Options:
		-id file   input_file for dominant model (SSGG annotation format, .psv for example)
		-ir file   input_file for recessive model (optional, if not defined, same file as for dominant model is assumed)
		-fi char   field_separator for input file (default = '|')
		-fo char   field_separator for output file (default = '\t')
		-s csv     comma-separated samples list (if not specified, all samples of input file are considered)
		-a csv     comma-separated list of affected_samples (if not specified, all samples of experiment are considered)
		-cd int    in the dominant mode, number of variants that a sample must at least contain to consider that it's affected in an heterozygous-compound mode (default value is 1, it must be >=1)
		-cr int    in the recessive mode, number of variants that a sample must at least contain to consider that it's affected in an heterozygous-compound mode (default value is 2, it must be >=2)
		-o string  output_folder (where global files and folder will be created. Default is cwd)
		-p string  prefix (included in all created file names)

		    Output files & folders:
		- 2 main folders are created, for the 2 disease inheritance modes (dominant or recessive). Inside these folders:
		    - sample_specific_results folder contains variants & genes lists for each individual sample and segregated into homozygous, heterozygous and mixed (compound) variants
		    - zygosity_specific_results folder contains variants & genes lists common to all affected samples segregated into homozozygous, heterozygous and compound (mixed genotypes over affected samples)
		- In the dominant folder, besides:
		    - dominant genes & files: variants either in homozygosis or heterozygosis for all affected samples
		- In the recessive folder, besides:
		    - recessive genes & files: 2 variants common to all affected samples, or a single variant in homozygosis common to all affected samples
		    - compound genes & files: at least \"c\" variants in each affected sample for each gene, not necessarily common to all affected samples
\n".localtime()."\n=================================================================\n\n";
}

if (not defined $input_recessive){
	$input_recessive = $input_dominant;
}
print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------

sub add_to_consequence_types_hash
# Function that adds a new variant object to a variants object hash, containing only variant objects that come from different consequence types of the same variant
# It checks if consequence type comes from the same variant, and if not all consequence types are first deleted and new one get stored
# A consequence type "number" is assigned, starting from 0 and incrementally
{
	my ($CTs_hash,$variant) = @_;
	# Store last variant and its different consequence types (to be recovered in case that hetero-compound mode determined)
	my $var_pos = 0;
	if (defined($$CTs_hash{0})){
		# Check last variants number:
		foreach (keys(%{$CTs_hash})){
			$var_pos ++;
		}
		# In case that a new variant is found (that is not a different consequence type of the same variant)
		if ($variant->compare($$CTs_hash{0})){
			foreach my $CT_num (keys(%{$CTs_hash})){
				delete($$CTs_hash{$CT_num});
			}
			$var_pos = 0;
		}
	}
	$$CTs_hash{$var_pos} = $variant;
}

sub main (){
	my $path = `pwd`;
	chomp $path;
	
	# Load input configuration
	my $config = analysis_info->get_config(psv_file => $input_dominant);
	
	# Get sample ids
	my @all_sample_ids = $config->get_sample_ids();
	my (@sample_ids,@affected_ids);
	if ($samples ne ""){
		@sample_ids = split(',',$samples);
		# Check if input samples are included in input file:
		foreach my $sample1 (@sample_ids){
			my $found = 0;
			foreach my $sample2 (@all_sample_ids){
				if ($sample1 eq $sample2){
					$found = 1;
				}
			}
			if (not $found){
				die ("ERROR: Sample $sample1 not included in input file $input_dominant!!");
			}
		}
	} else {
		@sample_ids = @all_sample_ids;
	}
	if ($affected eq ""){
		@affected_ids = @sample_ids;
	}else{
		@affected_ids = split(',',$affected);
	}
	# get the new configuration object:
	my $new_config = analysis_info->create_config(\@sample_ids,$config->{type},$sep_out);
	
	# IMPORTANT:
	# Since different filtering values may have been applied for both inheritance modes, script must run twice, one time with the dominant input file and other with the recessive input file (if different):
	foreach my $inheritance_mode (@inheritance_modes){
		# Open output files (hashes of file handles is defined)
		my (%output);
		mkdir($outfolder."/");
		mkdir($outfolder."/".$inheritance_mode);
		mkdir($outfolder."/".$inheritance_mode."/".$zygoresultspath);
		mkdir($outfolder."/".$inheritance_mode."/".$sampleresultspath);
		# Define headers of genes and variants files
		my %header;
		$header{'genes'} = "#HGNC_symbol".$sep_out."Total variants in gene from all affected samples".$sep_out."Variants in homozygosis in all affected samples".$sep_out."Variants in heterozygosis in all affected samples".$sep_out."Description\n";
		$header{'variants'} = $new_config->create_header()."\n";

		# In case of "|"-separated files, the extension of variants files is changed to .psv (for GeneSys compatibility)
		my $vars_ext = ".txt";
		if ($sep_out eq "|"){
			$vars_ext = ".psv";
		}
		# Create file handle for genes counts:
		my $gene_counter_file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/gene_counter.txt");
		$output{$inheritance_mode}{'gene_counter'} = $gene_counter_file_handle;
		#Gene counter header
		my $gene_counter_header = "#HGNC_symbol".$sep_out;			
		foreach my $sample (@affected_ids){
			#Gene counter header
			$gene_counter_header .= $sample.$sep_out;
		}
		chop $gene_counter_header;
		$gene_counter_header .= "\n";
		print {$output{$inheritance_mode}{'gene_counter'}} $gene_counter_header;
		
		foreach my $filetype (@file_types){
			my $file_handle;
			foreach my $inheritance (@inheritance_modes){
				if ($filetype eq "variants"){
					$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$prefix.$inheritance."_".$filetype.$vars_ext);
				} else {
					$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$prefix.$inheritance."_".$filetype.".txt");
				}
				$output{$inheritance}{$filetype} = $file_handle;
				print {$output{$inheritance}{$filetype}} $header{$filetype};
			}
			foreach my $zygos (@zygosity){
				if ($filetype eq "variants"){
					$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$zygoresultspath."/".$prefix.$filetype."_".$zygos.$vars_ext);
				} else {
					$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$zygoresultspath."/".$prefix.$filetype."_".$zygos.".txt");
				}
				$output{$zygoresultspath}{$filetype}{$zygos} = $file_handle;
				if ($filetype eq "genes"){
					print {$output{$zygoresultspath}{$filetype}{$zygos}} "#HGNC_symbol".$sep_out."Total $zygos variants in gene".$sep_out."Description\n";
				} else {
					print {$output{$zygoresultspath}{$filetype}{$zygos}} $header{'variants'};
				}
				foreach my $sample (@affected_ids){
					if ($filetype eq "variants"){
						$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$sampleresultspath."/".$sample."_".$filetype."_".$zygos.$vars_ext);
					} else {
						$file_handle = get_out_file_handle($outfolder."/".$inheritance_mode."/".$sampleresultspath."/".$sample."_".$filetype."_".$zygos.".txt");
					}
					$output{$sampleresultspath}{$sample}{$filetype}{$zygos} = $file_handle;
					if ($filetype eq "variants"){
						# Sample specific variants files only include one sample
						my @tmp_sample = ( $sample );
						my $sample_config = analysis_info->create_config(\@tmp_sample,$config->{type},$sep_out);
						my $tmp_header = $sample_config->create_header()."\n";
						print {$output{$sampleresultspath}{$sample}{$filetype}{$zygos}} $tmp_header;
					} elsif ($zygos eq "homozygous" or $zygos eq "heterozygous"){
						# Homozygous and heterozygous genes files must have a different header
						print {$output{$sampleresultspath}{$sample}{$filetype}{$zygos}} "#HGNC_symbol".$sep_out."Total $zygos variants in gene".$sep_out."Description\n";
					} else {
						print {$output{$sampleresultspath}{$sample}{$filetype}{$zygos}} $header{$filetype};
					}
				}
			}
		}
		
		# A hash containing important information is created, with following structure:
		# information{genes}{gene}{description}
		# information{genes}{gene}{samples}{sample}{homozygous}
		# information{genes}{gene}{samples}{sample}{heterozygous}
		# information{genes}{gene}{samples}{sample}{compound}
		# information{genes}{gene}{samples}{sample}{last_variant_CTs}
		# information{genes}{gene}{zygosity}{compound}
		# information{genes}{gene}{zygosity}{homozygous}
		# information{genes}{gene}{zygosity}{heterozygous}
		# information{genes}{gene}{all_samples_last_variant_CTs}
		# information{genes}{gene}{inheritance}{dominant}{all}
		# information{genes}{gene}{inheritance}{recessive}{all}
		my %information;
		
		# Open input file
		my $fh;
		if ($inheritance_mode eq "dominant"){
			$fh = get_in_file_handle($input_dominant);
		} else {
			$fh = get_in_file_handle($input_recessive);
		}
		
		# Parse each line of the annotation input file
		while (my $line = <$fh>){
			chomp($line);
			
			# Do not consider header or comment lines		
			if (not $line =~ m/^#/) {
				my $variant = analysis_info->get_fields_from_psv($line,$config);
				
				# Get gene name
				my $gene = $variant->GetGeneName();
				if (not defined($information{'genes'}{$gene})){
					$information{'genes'}{$gene} = {'description' => $variant->get_parameter("Gene_description"), 'samples' => {}};
				}
				
				# Variables for checking presence of variant in all samples depending on zygosity
				my $variant_in_all_samples = 1;
				my $homozygous_in_all_samples = 1;
				my $heterozygous_in_all_samples = 1;
				my $compound_in_all_samples = 1;
						
				foreach my $sample (@affected_ids){
					# Specific sample results can be printed out in this loop:
					if (genotype->IsVariant($variant->get_sample_parameter($sample,"Genotype"))){
						
						# Create genes tags (if not yet created)
						if (not defined($information{'genes'}{$gene}{'samples'}{$sample})){
							$information{'genes'}{$gene}{'samples'}{$sample} = {'homozygous' => 0,'heterozygous' => 0,'compound' => 0};
						}
						
						# Check if new variant position has changed from last variant for this sample
						my $position_changed = 0;
						if (not defined($information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'})){
							$position_changed = 1;
							$information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'} = {};
						} elsif ($variant->compare($information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'}{0})){
							$position_changed = 1;
						}
						
						# Config for sample-specific file: 
						my @tmp_sample = ( $sample );
						my $sample_config = analysis_info->create_config(\@tmp_sample,$config->{type},$sep_out);
						
						# Update per sample gene counts and print to sample specific files
						my $var_zygosity = genotype->Zygosity($variant->get_sample_parameter($sample,"Genotype"));
											
						if ($var_zygosity == 1){
							if ($position_changed){
								$information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'} ++;
							}
							$homozygous_in_all_samples = 0;
							print {$output{$sampleresultspath}{$sample}{'variants'}{'heterozygous'}} $variant->put_fields_into_psv($sample_config);
						} elsif ($var_zygosity == 2){
							if ($position_changed){
								$information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'} ++;
							}
							$heterozygous_in_all_samples = 0;
							print {$output{$sampleresultspath}{$sample}{'variants'}{'homozygous'}} $variant->put_fields_into_psv($sample_config);
						}
						
						# If at least 2 variants, then variant is also stored in 'compound' file
						my $gene_variants = $information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'} + $information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'};
						if ($gene_variants > 1){
							if ($position_changed){
								if ($gene_variants == 2){
									$information{'genes'}{$gene}{'samples'}{$sample}{'compound'} ++;
									foreach my $ct (keys(%{$information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'}})){
										print {$output{$sampleresultspath}{$sample}{'variants'}{'compound'}} $information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'}{$ct}->put_fields_into_psv($sample_config);
									}
								}
								$information{'genes'}{$gene}{'samples'}{$sample}{'compound'} ++;
							}
							print {$output{$sampleresultspath}{$sample}{'variants'}{'compound'}} $variant->put_fields_into_psv($sample_config);
						} else {
							$compound_in_all_samples = 0;
						}
											
						# Store last variant (to be recovered in case that compound mode determined)
						add_to_consequence_types_hash($information{'genes'}{$gene}{'samples'}{$sample}{'last_variant_CTs'},$variant);
					} else {
						$variant_in_all_samples = 0;
					}
				}
				
				# When variant is found in all affected samples, gene is considered as candidate
				if ($variant_in_all_samples){
					# Check if new variant position has changed from last variant for this sample
					my $position_changed = 0;
					if (not defined($information{genes}{$gene}{'all_samples_last_variant_CTs'})){
						$position_changed = 1;
						$information{genes}{$gene}{'all_samples_last_variant_CTs'} = {};
					} elsif ($variant->compare($information{genes}{$gene}{'all_samples_last_variant_CTs'}{0})){
						$position_changed = 1;
					}
					# Create genes tags (if not yet created)
					if (not defined($information{genes}{$gene}{'zygosity'})){
						$information{genes}{$gene}{'zygosity'} = {'homozygous'=> 0,
							'heterozygous' => 0,
							'compound' => 0};
					}
					if (not defined($information{genes}{$gene}{'inheritance'})){
						$information{genes}{$gene}{'inheritance'} = {'dominant' => {'all' => 0},
							'recessive' => {'all' => 0}};
					}
					
					# Count and print variants for zygosity and inheritance output files
					my $gene_variants = $information{'genes'}{$gene}{'inheritance'}{'dominant'}{'all'};
					my $gene_rec_variants = $information{'genes'}{$gene}{'inheritance'}{'recessive'}{'all'};
					
					# Print output to variants files
					print {$output{'dominant'}{'variants'}} $variant->put_fields_into_psv($new_config);
					# for the recessive model, may be there was already an heterozygous variant not added (one hetero variant is not enough for recessive model)
					if ($position_changed){
						$information{'genes'}{$gene}{'inheritance'}{'dominant'}{'all'} ++;
						# Store and count in corresponding files
						if ($homozygous_in_all_samples){
							$information{'genes'}{$gene}{'zygosity'}{'homozygous'} ++;
							print {$output{$zygoresultspath}{'variants'}{'homozygous'}} $variant->put_fields_into_psv($new_config);
						} elsif ($heterozygous_in_all_samples){
							$information{'genes'}{$gene}{'zygosity'}{'heterozygous'} ++;
							print {$output{$zygoresultspath}{'variants'}{'heterozygous'}} $variant->put_fields_into_psv($new_config);
						} else {
							$information{'genes'}{$gene}{'zygosity'}{'compound'} ++;
							print {$output{$zygoresultspath}{'variants'}{'compound'}} $variant->put_fields_into_psv($new_config);
						}
						if ($gene_variants > 0 or ($gene_variants == 0 and $homozygous_in_all_samples == 1)){
							if ($gene_variants == 1 and $gene_rec_variants == 0){
								$information{'genes'}{$gene}{'inheritance'}{'recessive'}{'all'} ++;
								foreach my $ct (keys(%{$information{'genes'}{$gene}{'all_samples_last_variant_CTs'}})){
									print {$output{'recessive'}{'variants'}} $information{'genes'}{$gene}{'all_samples_last_variant_CTs'}{$ct}->put_fields_into_psv($new_config);
								}
							}
							$information{'genes'}{$gene}{'inheritance'}{'recessive'}{'all'} ++;
							print {$output{'recessive'}{'variants'}} $variant->put_fields_into_psv($new_config);
						}
					} else {
						if ($homozygous_in_all_samples){
							print {$output{$zygoresultspath}{'variants'}{'homozygous'}} $variant->put_fields_into_psv($new_config);
						} elsif ($heterozygous_in_all_samples){
							print {$output{$zygoresultspath}{'variants'}{'heterozygous'}} $variant->put_fields_into_psv($new_config);
						}
						if ($gene_variants > 1 or ($gene_variants == 1 and $homozygous_in_all_samples == 1)){
							print {$output{'recessive'}{'variants'}} $variant->put_fields_into_psv($new_config);
						}
					}
					# Store last variant and its different consequence types (to be recovered in case that hetero-compound mode determined)
					add_to_consequence_types_hash($information{'genes'}{$gene}{'all_samples_last_variant_CTs'},$variant);
				}
			}
		}
		
		# close input file
		$fh->close();
		
		# Print information to genes files
		foreach my $gene (sort(keys(%{$information{'genes'}}))){
			foreach my $inheritance (@inheritance_modes){
				my $variants = $information{'genes'}{$gene}{'inheritance'}{$inheritance}{'all'};
				if (defined($variants) and $variants > 0){
					my $string = $gene.$sep_out.$variants.$sep_out.$information{'genes'}{$gene}{'zygosity'}{'homozygous'}.$sep_out.$information{'genes'}{$gene}{'zygosity'}{'heterozygous'}.
						$sep_out.$information{'genes'}{$gene}{'description'}."\n";
					print {$output{$inheritance}{'genes'}} $string;
				}
			}
			
			#Gene counter
			my $gene_counter_string = "$gene";
			my $variant_in_any_sample = 0;
			foreach my $sample (@affected_ids){
				my ($homozygous_sum, $heterozygous_sum, $affected_sum);			
				if (not defined($information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'})){
					$information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'} = 0;
				}
				$homozygous_sum = $information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'};
				
				if (not defined($information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'})){
					$information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'} = 0;
				}
				if (not defined($information{'genes'}{$gene}{'samples'}{$sample}{'compound'})){
					$information{'genes'}{$gene}{'samples'}{$sample}{'compound'} = 0;
				}
				$heterozygous_sum = $information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
				$affected_sum = $homozygous_sum + $heterozygous_sum;
				if ($affected_sum > 0){
					$variant_in_any_sample = 1;
				}
				$gene_counter_string .= $sep_out.$affected_sum;
			}
			if ($variant_in_any_sample){
				print {$output{$inheritance_mode}{'gene_counter'}} $gene_counter_string."\n";
			}
			
			foreach my $zygos (@zygosity){
				my $string = "";
				if (defined($information{'genes'}{$gene}{'zygosity'}{$zygos}) and $information{'genes'}{$gene}{'zygosity'}{$zygos} > 0){
					$string = $gene.$sep_out.$information{'genes'}{$gene}{'zygosity'}{$zygos};
				}
				if ($string ne ""){
					$string .= $sep_out.$information{'genes'}{$gene}{'description'}."\n";
					print {$output{$zygoresultspath}{'genes'}{$zygos}} $string;
				}
				foreach my $sample (@affected_ids){
					$string = "";
					if (defined($information{'genes'}{$gene}{'samples'}{$sample})){
						my $homozygous_vars = 0;
						if (defined($information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'})){
							$homozygous_vars = $information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'};
						}
						my $heterozygous_vars = 0;
						if (defined($information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'})){
							$heterozygous_vars = $information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
						}
						if ($zygos eq "homozygous"){
							if ($information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'} > 0){
								$string = $gene.$sep_out.$information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'};
							}
						} elsif ($zygos eq "heterozygous"){
							if ($information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'} > 0){
								$string = $gene.$sep_out.$information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
							}
						} elsif ($zygos eq "compound"){
							my $all_vars = $information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'} + $information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
							if ($all_vars > 1){
								$string = $gene.$sep_out.$all_vars.$sep_out.$gene.$sep_out.$information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'}.$sep_out.$information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
							}
						}
						if ($string ne ""){
							$string .= $sep_out.$information{'genes'}{$gene}{'description'}."\n";
							print {$output{$sampleresultspath}{$sample}{'genes'}{$zygos}} $string;
						}
					}
				}
			}
		}
		
		# TODO: now all steps are performed twice, even file output. At the end, either the dominant or recessive file must be deleted
		
		# Close files
		foreach my $filetype (@file_types){
			foreach my $inheritance (@inheritance_modes){
				$output{$inheritance}{$filetype}->close();
				if ($inheritance_mode ne $inheritance){
					if ($filetype eq "variants"){
						unlink($outfolder."/".$inheritance_mode."/".$prefix.$inheritance."_".$filetype.$vars_ext);
					} else {
						unlink($outfolder."/".$inheritance_mode."/".$prefix.$inheritance."_".$filetype.".txt");
					}
				}
			}
			foreach my $zygos (@zygosity){
				$output{$zygoresultspath}{$filetype}{$zygos}->close();
				foreach my $sample (@affected_ids){
					$output{$sampleresultspath}{$sample}{$filetype}{$zygos}->close();
				}
			}
		}
		
		# Compound modes (where genes with variants in different positions are also included as candidates) must also be considered
		# --------------

		# Open again the input file and parse it:
		my $compound;
		if ($inheritance_mode eq "dominant"){
			$fh = get_in_file_handle($input_dominant);
			$compound = $compound_dominant;
		} else {
			$fh = get_in_file_handle($input_recessive);
			$compound = $compound_recessive;
		}
		print "INFO: Inheritance mode $inheritance_mode, minimum number of compound-variants per gene is $compound\n";
		
		# Open output files and print headers
		my $compound_vars_fh = get_out_file_handle($outfolder."/".$inheritance_mode."/".$prefix."compound_variants.txt");
		print $compound_vars_fh $header{'variants'};
		my $compound_genes_fh = get_out_file_handle($outfolder."/".$inheritance_mode."/".$prefix."compound_genes.txt");
		print $compound_genes_fh "#HGNC_symbol".$sep_out."Total variants in gene from all affected samples".$sep_out."Variants present in all affected samples".$sep_out."Variants in homozygosis in all affected samples".$sep_out."Variants in heterozygosis in all affected samples".$sep_out."Description\n";
			
		# Parse input file
		my $last_variant;
		while (my $line = <$fh>){
			chomp($line);
			# Do not consider header or comment lines		
			if (not $line =~ m/^#/) {
				my $variant = analysis_info->get_fields_from_psv($line,$config);
				# Get gene name
				my $gene = $variant->GetGeneName();
				
				# check if all affected samples have at least "$compound" variants in the variants gene, and if so, output to compound file
				my $compound_var = 1;
				my $is_variant_in_any_sample = 0;
				foreach my $sample (@affected_ids){
					# Obtain number of variants per gene in the selected sample
					my $homozygous_vars = 0;
					if (defined($information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'})){
						$homozygous_vars = $information{'genes'}{$gene}{'samples'}{$sample}{'homozygous'};
					}
					my $heterozygous_vars = 0;
					if (defined($information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'})){
						$heterozygous_vars = $information{'genes'}{$gene}{'samples'}{$sample}{'heterozygous'};
					}
					my $sample_gene_variant_alleles = $heterozygous_vars + 2 * $homozygous_vars;
					if (not defined($sample_gene_variant_alleles)){
						$sample_gene_variant_alleles = 0;
					}
					if ($sample_gene_variant_alleles < $compound){
						$compound_var = 0;
					}
					if (genotype->IsVariant($variant->get_sample_parameter($sample,"Genotype"))){
						$is_variant_in_any_sample = 1;
					}
				}
				if ($compound_var and $is_variant_in_any_sample){
					if (not defined($information{'genes'}{$gene}{'inheritance'}{'recessive'}{'present_at_any_affected_sample'})){
						$information{'genes'}{$gene}{'inheritance'}{'recessive'}{'present_at_any_affected_sample'} = 1;
						$last_variant = $variant;
					} elsif ($variant->compare($last_variant)) {
						$information{'genes'}{$gene}{'inheritance'}{'recessive'}{'present_at_any_affected_sample'} ++;
						$last_variant = $variant;
					}
					print $compound_vars_fh $variant->put_fields_into_psv($new_config);
				}
			}
		}
		
		foreach my $gene (sort(keys(%{$information{'genes'}}))){
			# Get variants present in all samples
			my $vars_present_all_samps = 0;
			foreach my $zygos (@zygosity){
				if (defined($information{'genes'}{$gene}{'zygosity'}{$zygos})){
					$vars_present_all_samps += $information{'genes'}{$gene}{'zygosity'}{$zygos};
				} elsif (defined($information{'genes'}{$gene}{'zygosity'})) {
					$information{'genes'}{$gene}{'zygosity'}{$zygos} = 0;
				} else {
					$information{'genes'}{$gene}{'zygosity'} = { $zygos => 0};
				}
			}
			if (defined($information{'genes'}{$gene}{'inheritance'}{'recessive'}{'present_at_any_affected_sample'})){
				print $compound_genes_fh $gene.$sep_out.$information{'genes'}{$gene}{'inheritance'}{'recessive'}{'present_at_any_affected_sample'}.$sep_out.$vars_present_all_samps.$sep_out.$information{'genes'}{$gene}{'zygosity'}{'homozygous'}.$sep_out.$information{'genes'}{$gene}{'zygosity'}{'heterozygous'}.$sep_out.$information{'genes'}{$gene}{'description'}."\n";
			}
		}
		
		# Close files
		$compound_vars_fh->close();
		$compound_genes_fh->close();
		$fh->close();
	}
}

exit;