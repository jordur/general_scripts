#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa, Arbol
# @Contributors: 
########################################################

package var_calling;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use NGS::essay; 
use file_handle;
use submit;
use Exporter;
use variants_info;
use genomics::genotype;
use File::Spec;

# -----------------------
# -- Global variables  --
# -----------------------
my $module = "variant_calling";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub variant_calling {
	my ($class, $essay, $id) = @_;
	print "Launching variant calling jobs for essay $essay->{name}... ".localtime()."\n";
	
	# Illumina analyses are performed only with PE reads, while SOLiD add also F3
	my @types;
	if ($essay->is_solid() eq 'true'){
		@types = ("F3", "PE");
	} else {
		@types = ("PE");
	}
	
	my $var_id;
	
	# The _sample_replicates.txt file is needed for passing the sample names argument to VarScan. Since many VarScan variant calling processes 
	# can be performed simultaneously, this file is better created here
	my $sample_replicates_list = $essay->add_file("trash","variants",$essay->{'name'}."_sample_replicates.txt");
	# For avoiding problems in headers when collecting, NO type (PE,F3) information will be added to sample names:
	my $replicates = join("\t",@{$essay->get_replicates()});
	$essay->create_tree();
	my $string = `echo \"$replicates\" > $sample_replicates_list\n`;
	
	foreach my $type (@types){
		my $var_call_id;
		
		# If target sequencing panel is big (exomes), a chromosome-split of the mapping bam is recommended
		if (defined($essay->{target_reference}{chr_split}) and $essay->{target_reference}{chr_split} eq "yes"){

			# First launch all split jobs
			my $chr_split_id;
			foreach my $sample (@{$essay->{samples_order}}){
				foreach my $replicate (sort(keys(%{$essay->{samples}{$sample}{replicates}}))){
					my @bams = @{$essay->get_bams(sample=>$sample,replicate=>$replicate,level=>"replicate")};
					my $split_id = &split_chrs($essay,$id,$bams[0],$sample,$replicate,$type);
					$chr_split_id = add_job_ids($chr_split_id,$split_id);
				}
			}
			
			# Once ready with splitting, variant calling can be performed
			my @chroms = @{genotype->Chromosomes()};
			foreach my $chr (@chroms){
				my $gatk_id = &gatk ($essay,$chr_split_id,$type,$chr);
				my $varsc_id = &varscan ($essay,$chr_split_id,$type,$chr);
						
				# Next job must wait for all callers
				$var_call_id = add_job_ids($var_call_id,$gatk_id);
				$var_call_id = add_job_ids($var_call_id,$varsc_id);
			}
			
			# Now join all chromosomes
			my $chr_join_gatk_id = &join_chrs($essay,$var_call_id,"gatk",$type);
			my $chr_join_varscan_id = &join_chrs($essay,$var_call_id,"varscan",$type);
			
			# Next job must wait for all callers
			$var_call_id = add_job_ids($var_call_id,$chr_join_gatk_id);
			$var_call_id = add_job_ids($var_call_id,$chr_join_varscan_id);
		} else {
			my $gatk_id = &gatk ($essay,$id,$type);
			my $varsc_id = &varscan ($essay,$id,$type);

			# Next job must wait for all callers
			$var_call_id = add_job_ids($var_call_id,$gatk_id);
			$var_call_id = add_job_ids($var_call_id,$varsc_id);
		}
	
		# At the moment bioscope small indels calling is performed separately and orchestrated by BioScope
		# These operations must be performed before launching MeLISA, and no job_id is therefore returned
		if ($essay->is_solid() eq 'true'){
			&bioscope_gff ($essay,$var_call_id,$type);
		}
		$var_id = add_job_ids($var_id,$var_call_id); 
	}
	
	return &collecting ($essay,$var_id);
}

sub bioscope_gff {
	my ($essay,$previous_jobs,$type) = @_;
	my $vcfs = "";
	my $output = $essay->add_file("trash","variants",$essay->{name}."_gff_".$type);
	foreach my $sample (keys(%{$essay->get_parameter("samples")})){
		foreach my $replicate (keys(%{$essay->get_parameter(("samples",$sample,"replicates"))})){
			my $bam;
			if ($type eq "F3"){
				$bam = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge_F3.bam");
			} else {
				$bam = $essay->add_file("trash","mapping",$sample,$replicate,$sample."_".$replicate."_merge.bam");
			}
			my $output_bioscope = $essay->mkdir("trash","variants",$sample,$replicate,"bioscope_small_indels_".$type);
			my $vcf = $essay->add_file("trash","variants",$sample,$replicate,"indels_".$type."_bioscope.vcf");
			my $output_path = $essay->mkdir("trash","variants",$sample,$replicate);
			my $vcf_in_range = "indels_in_range_".$type."_bioscope.vcf";
			my $gff = $essay->add_file("trash","variants",$sample,$replicate,"bioscope_small_indels_".$type,"output","smallindel","New_Experiment.gff");
			my $target = $essay->{target_reference}{target};
			my $job_file = $essay->add_file("jobs","variants",$sample,$replicate,"bioscope_small_indels_".$type);
			my $job_file2 = $essay->add_file("jobs","variants",$sample,$replicate,"vcf_bioscope_small_indels_".$type);
			my $log_dir = $essay->mkdir("logs","variants",$sample,$replicate);
			
			# Create involved folders
			$essay->create_tree();
			
			my $string = "sg_bioscope_analyzer.pl -pipe smallindels -seq $type -bam $bam -outdir $output_bioscope -run";
			my $job = launch_task($job_file,'',"",$string,$log_dir,$essay,$module);
			
			$string = "sg_parse2JM_format.sh -s $sample.$replicate -o $vcf $gff
			sg_extraer_indels_en_rango.sh $target $vcf $vcf_in_range $output_path";
			$job = launch_essay_job ($job_file2, "", $previous_jobs, $string, $log_dir, $essay, $module);
			
			$vcfs = $vcfs.$output_path."/".$vcf_in_range.",";
		}
	}
	$vcfs =~ s/,$//;
	my $job_file = $essay->add_file("jobs","variants","bioscope_small_indels_".$type);
	my $log_dir = $essay->mkdir("logs","variants");
	my $string = "sg_parsing_indels.pl -v $vcfs -o $output";
	my $job = launch_essay_job ($job_file, "", $previous_jobs, $string, $log_dir, $essay, $module);
}

sub split_chrs {
	# This function is currently a mapping function
	
	my ($essay,$previous_jobs,$bam,$sample,$replicate,$type) = @_;
	
	# Define involved files
	my $job_file = $essay->add_file("jobs","variants",$sample,$replicate,"split_".$type);
	my $chr_split;
	if (($essay->is_solid() eq 'true') and ($type eq "F3")){
		$chr_split = $essay->add_file("trash","variants",$sample,$replicate,$sample."_".$replicate."_split".$type);
	} else {
		$chr_split = $essay->add_file("trash","variants",$sample,$replicate,$sample."_".$replicate."_split");
	}
	my $log_dir = $essay->mkdir("logs","variants",$sample,$replicate);
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "bamtools split -in $bam -reference -stub $chr_split
	for i in $chr_split*.bam; do samtools index \$i; done";
	my $job = launch_essay_job ($job_file, "", $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub join_chrs
# using variants_info class, perform the joining of the chromosomes files
{
	my ($essay,$previous_jobs,$var_caller,$type) = @_;
	
	my @chroms;
	my $dir = $essay->get_parameter("target_reference","target_chrs");
	opendir(DIR, $dir) or die $!;
    while (my $file = readdir(DIR)) {
    	my $path = File::Spec->catfile($dir, $file);
    	# Use a regular expression to ignore files beginning with a period
        next if ($file =~ m/^\./);
        next if (-z $path);
		$file =~ s/chr//;
		$file =~ s/.bed//;
		push(@chroms,$file);
    }
    closedir(DIR);
    
	my $analysis_name = $essay->{'target_reference'}{'name'};
	#Modify chr to only use 13 and 17 if we're on Multiplicom BRCAs
	if ($analysis_name eq 'brcas_multiplicom'){
		@chroms = (13, 17);
	}
	
	my $files = "";
	foreach my $chr (@chroms){
		my $file = $essay->add_file("trash","variants","split_".$type.".REF_chr".$chr."_".$var_caller."_normalized.vcf");
		$files .= "--variant:".$var_caller."_chr".$chr." ".$file." ";
	}
	
	my $reference = $essay->{references}{fasta};  # Whole genome reference, otherwise variant callers don't work properly
	
	my $job_file =  $essay->add_file("jobs","variants","join_chrs_".$type."_".$var_caller);
	my $output;
	if (($essay->is_solid() eq 'true') and ($type eq "F3")){
		$output = $essay->add_file("trash","variants",$essay->{name}."_".$var_caller."_".$type."_normalized.vcf");
	} else {
		$output = $essay->add_file("trash","variants",$essay->{name}."_".$var_caller."_normalized.vcf");
	}
	my $log_dir = $essay->mkdir("logs","variants");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "gatk_toolkit -R $reference -T CombineVariants $files -o $output -mergeInfoWithMaxAC -printComplexMerges -assumeIdenticalSamples";
	my $job = launch_essay_job ($job_file, "", $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub gatk {
	my ($essay,$previous_jobs,$type,$chr) = @_;

	# Define involved files
	my $job_file;
	my $reference = $essay->{'references'}{'fasta'}; # Whole genome reference, otherwise variant callers don't work properly
	my $analysis_name = "";
	if (defined($essay->{'target_reference'}{'name'})){
		$analysis_name = $essay->{'target_reference'}{'name'};
	}
	my ($bams_string,$gatk,$gatk_norm,$gatk_raw_variants);
	my $target = "";
	my $options = "";
	my $string = "";
	my $dcov = '10000';
	my $threads = 1;
	
	if (($essay->is_solid() eq 'true') and ($type eq "F3")){
		if (defined($chr)){
			my @bams = @{$essay->get_bams(type=>$type,level=>"replicate",chr=>$chr)};
			foreach my $bam (@bams){
				$bams_string = $bams_string . " -I $bam";
			}
			$job_file =  $essay->add_file("jobs","variants","gatk_var_calling_F3_chr".$chr);
			$gatk_raw_variants = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_gatk_raw");
			$gatk = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_gatk.vcf");
			$gatk_norm = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_gatk_normalized.vcf");
		} else {
			my @bams = @{$essay->get_bams(type=>$type,level=>"replicate")};
			foreach my $bam (@bams){
				$bams_string = $bams_string . " -I $bam";
			}
			$job_file =  $essay->add_file("jobs","variants","gatk_var_calling_F3");
			$gatk_raw_variants = $essay->add_file("trash","variants",$essay->{name}."_gatk_raw_F3");
			$gatk = $essay->add_file("trash","variants",$essay->{name}."_gatk_F3.vcf");
			$gatk_norm = $essay->add_file("trash","variants",$essay->{name}."_gatk_F3_normalized.vcf");
			$options = "-pe smp 4";
			$threads = 4;
		}
	} else {
		if (defined($chr)){
			my @bams = @{$essay->get_bams(type=>$type,level=>"replicate",chr=>$chr)};
			foreach my $bam (@bams){
				$bams_string = $bams_string . " -I $bam";
			}
			$job_file =  $essay->add_file("jobs","variants","gatk_var_calling_PE_chr".$chr);
			$gatk_raw_variants = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_gatk_raw");
			$gatk = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_gatk.vcf");
			$gatk_norm = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_gatk_normalized.vcf");
		} else {
			my @bams = @{$essay->get_bams(type=>$type,level=>"replicate")};
			foreach my $bam (@bams){
				$bams_string = $bams_string . " -I $bam";
			}
			$job_file =  $essay->add_file("jobs","variants","gatk_var_calling");
			$gatk_raw_variants = $essay->add_file("trash","variants",$essay->{name}."_gatk_raw");
			$gatk = $essay->add_file("trash","variants",$essay->{name}."_gatk.vcf");
			$gatk_norm = $essay->add_file("trash","variants",$essay->{name}."_gatk_normalized.vcf");
			$options = "-pe smp 4";
			$threads = 4;
		}
	}
	if (defined($chr)){
		if (defined($essay->{'target_reference'}{'target_chrs'})){
			$target = "-L ".$essay->{'target_reference'}{'target_chrs'}."/chr".$chr.".bed";
		}
	} else {
		if (defined($essay->{'target_reference'}{'target'})){
			$target = "-L ".$essay->{'target_reference'}{'target'};
		}
	}
	# For avoiding problems in headers when collecting, NO type (PE,F3) information will be added to sample names:
	my $replicates = join(',',@{$essay->get_replicates()});
	
	my $log_dir = $essay->mkdir("logs","variants");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	# TODO: Remove the -rf BadCigar as soon as possible
	# TODO: Add dcov as a JSON input parameter, modify matraz config and script.
	#If we detect brcas_multiplicom we change dcov to launch GATK
	if ($analysis_name eq 'brcas_multiplicom'){
		$dcov = '10000';
		$options = "-pe smp 8";
		$threads = 8;
	}
	$string = "gatk_toolkit -T UnifiedGenotyper -rf BadCigar -R $reference -nt 1 -nct $threads $bams_string -o $gatk -glm BOTH -dcov $dcov $target --min_base_quality_score 0 -minIndelCnt 3 -stand_call_conf 0 -stand_emit_conf 0 -minIndelFrac 0\n";
	$string .= "normalize_variants_vcf.py -i $gatk -o $gatk_norm\n";
	my $job = launch_essay_job ($job_file, $options, $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub varscan {
	my ($essay,$previous_jobs,$type,$chr) = @_;

	my $max_depth_per_bam = '';
	# Define involved files
	my $job_file;
	my $reference = $essay->{references}{fasta}; # Whole genome reference, otherwise variant callers don't work properly
	my $analysis_name = "";
	if (defined($essay->{target_reference}{name})){
		$analysis_name = $essay->{target_reference}{name};
	}
	my ($varscan,$varscan_raw_variants,$samparam,$mpileup,$cns,$varscan_norm);
	my $target = "";
	if (($essay->is_solid() eq 'true') and ($type eq "F3")){
		if (defined($chr)){
			$job_file =  $essay->add_file("jobs","variants","varscan_var_calling_F3_chr".$chr);
			$mpileup = $essay->add_file("trash","variants","varscan_mpileup_F3_chr".$chr);
			$cns = $essay->add_file("trash","variants","varscan_cns_F3_chr".$chr);
			$varscan_raw_variants = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_varscan_raw");
			$varscan = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_varscan.vcf");
			$varscan_norm = $essay->add_file("trash","variants","split_F3.REF_chr".$chr."_varscan_normalized.vcf");
		} else {
			$job_file =  $essay->add_file("jobs","variants","varscan_var_calling_F3");
			$mpileup = $essay->add_file("trash","variants",$essay->{name}."_varscan_mpileup_F3");
			$cns = $essay->add_file("trash","variants",$essay->{name}."_varscan_cns_F3");
			$varscan_raw_variants = $essay->add_file("trash","variants",$essay->{name}."_varscan_raw_F3");
			$varscan = $essay->add_file("trash","variants",$essay->{name}."_varscan_F3.vcf");
			$varscan_norm = $essay->add_file("trash","variants",$essay->{name}."_varscan_F3_normalized.vcf");
		}
	} else {
		if (defined($chr)){
			$job_file =  $essay->add_file("jobs","variants","varscan_var_calling_PE_chr".$chr);
			$mpileup = $essay->add_file("trash","variants","varscan_mpileup_PE_chr".$chr);
			$cns = $essay->add_file("trash","variants","varscan_cns_PE_chr".$chr);
			$varscan_raw_variants = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_varscan_raw");
			$varscan = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_varscan.vcf");
			$varscan_norm = $essay->add_file("trash","variants","split_PE.REF_chr".$chr."_varscan_normalized.vcf");
		} else {
			$job_file =  $essay->add_file("jobs","variants","varscan_var_calling");
			$mpileup = $essay->add_file("trash","variants",$essay->{name}."_varscan_mpileup");
			$cns = $essay->add_file("trash","variants",$essay->{name}."_varscan_cns");
			$varscan_raw_variants = $essay->add_file("trash","variants",$essay->{name}."_varscan_raw");
			$varscan = $essay->add_file("trash","variants",$essay->{name}."_varscan.vcf");
			$varscan_norm = $essay->add_file("trash","variants",$essay->{name}."_varscan_normalized.vcf");
		}
	}
	if (defined($chr)){
		if (defined($essay->{target_reference}{target_chrs})){
			$target = "-l ".$essay->{target_reference}{target_chrs}."/chr".$chr.".bed";
		}
		my @bams = @{$essay->get_bams(type=>$type,level=>"replicate",chr=>$chr)};
		foreach my $bam (@bams){
			$samparam = $samparam . " $bam";
		}		
	} else {
		if (defined($essay->{target_reference}{target})){
			$target = "-l ".$essay->{target_reference}{target};
		}
		my @bams = @{$essay->get_bams(type=>$type,level=>"replicate")};
		foreach my $bam (@bams){
			$samparam = $samparam . " $bam";
		}
	}
	
	my $sample_replicates_list = $essay->add_file("trash","variants",$essay->{'name'}."_sample_replicates.txt");
	my $log_dir = $essay->mkdir("logs","variants");

	# Create involved folders
	$essay->create_tree();
	
	#If we detect brcas_multiplicom we add max depth per bam option
	if ($analysis_name eq 'brcas_multiplicom'){
		$max_depth_per_bam = '-d 100000';
	}
	
	# Launch jobs
	my $string = "samtools mpileup $max_depth_per_bam -A -D -S -m 1 -F 0.001 -P Illumina $target -f $reference $samparam > $mpileup\n";
	$string .= "count_lines=\`wc -l $mpileup | awk \'{print \$1}\'\`\n";
	$string .= "if [[ \$count_lines > 0 ]]; then\n";
	# VarScan was apparently creating a good amount of noise in variant calls. Therefore a min-coverage of 20 and min-var-freq of 0.1 are set
	$string .= "\tVarScan mpileup2cns $mpileup --min-avg-qual 0 --min-coverage 3 --min-freq-for-hom 0.85 --min-var-freq 0.1 --p-value 0.99 --strand-filter 0 --variants --output-vcf 1 --vcf-sample-list $sample_replicates_list > $cns\n";
	$string .= "\trepairing_VarScan_vcf.py -i $cns -o $varscan\n";
	$string .= "\tnormalize_variants_vcf.py -i $varscan -o $varscan_norm\n";
	$string .= "fi";
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub collecting {
	my ($essay,$previous_jobs) = @_;

	# Define involved files
	my $job_file = $essay->add_file("jobs","variants","collecting");
	my $collect = $essay->add_file("trash","variants",$essay->{name}."_collect.vcf");
	my $collect_snps = $essay->add_file("trash","variants",$essay->{name}."_collect_snps.vcf");
	my $collect_snps_filtered = $essay->add_file("trash","variants",$essay->{name}."_collect_snps_filtered.vcf");
	my $collect_indels = $essay->add_file("trash","variants",$essay->{name}."_collect_indels.vcf");
	my $collect_indels_filtered = $essay->add_file("trash","variants",$essay->{name}."_collect_indels_filtered.vcf");
	my $final = $essay->add_file("analysis","variants",$essay->{name}."_collect.vcf");
	my $gatk = $essay->add_file("trash","variants",$essay->{name}."_gatk_normalized").".vcf";
	my $varscan = $essay->add_file("trash","variants",$essay->{name}."_varscan_normalized").".vcf";
	my $snp_af_threshold = genotype->get_af_threshold("SNP");
	my $indel_af_threshold = genotype->get_af_threshold("INDEL");
	my $depth_threshold = genotype->get_depth_threshold();
	my $reference = $essay->{references}{fasta};  # Whole genome reference, otherwise variant callers don't work properly
	my $target = "";
	if (defined($essay->{target_reference}{target})){
		$target = "-L ".$essay->{target_reference}{target};
	}
	
	# BioScope alternative is deprecated and shouldn't work...
	#my $bp = "";
	#if (defined($essay->get_path("trash","variants",$essay->{name}."_gff_PE"))){
	#	print "INFO: Including variants coming from BioScope PE (".$essay->get_path("trash","variants",$essay->{name}."_gff_PE").")\n";
	#	$bp = "-bp ".$essay->get_path("trash","variants",$essay->{name}."_gff_PE").".snvs";
	#}
	#my $gf = "";
	#if (defined($essay->get_path("trash","variants",$essay->{name}."_gatk_F3"))){
	#	print "INFO: Including variants coming from Gatk F3 (".$essay->get_path("trash","variants",$essay->{name}."_gatk_F3").".snvs)\n";
	#	$gf = "-gf ".$essay->get_path("trash","variants",$essay->{name}."_gatk_F3").".snvs";
	#}
	#my $vf = "";
	#if (defined($essay->get_path("trash","variants",$essay->{name}."_varscan_F3"))){
	#	print "INFO: Including variants coming from VarScan F3 (".$essay->get_path("trash","variants",$essay->{name}."_varscan_F3").".snvs)\n";
	#	$vf = "-vf ".$essay->get_path("trash","variants",$essay->{name}."_varscan_F3").".snvs";
	#}
	#my $bf = "";
	#if (defined($essay->get_path("trash","variants",$essay->{name}."_gff_F3"))){
	#	print "INFO: Including variants coming from BioScope F3 (".$essay->get_path("trash","variants",$essay->{name}."_gff_F3").")\n";
	#	$bf = "-bf ".$essay->get_path("trash","variants",$essay->{name}."_gff_F3").".snvs";
	#}
	
	my $log_dir = $essay->mkdir("logs","variants");
	
	# Create involved folders
	$essay->create_tree();
	
	# Launch jobs
	my $string = "collecting_vcf_variants.py -i $gatk -i $varscan -o $collect\n";
	if ($essay->get_parameter("modules","variant_calling","parameters","filter") eq "yes"){
		# Temporary changes. Snps and Indels directly splitted out from collect file in collecting_vcf_variants script
		$string .= "mv ${collect}_snps.vcf $collect_snps\n";
		$string .= "mv ${collect}_indels.vcf $collect_indels\n";
		#$string .= "gatk_toolkit -T SelectVariants -R $reference -V $collect $target -selectType SNP -o $collect_snps\n";
		#$string .= "gatk_toolkit -T SelectVariants -R $reference -V $collect $target -selectType INDEL -o $collect_indels\n";
		$string .= "gatk_toolkit -T VariantFiltration -R $reference -V $collect_snps --filterExpression \"SFREQ \< $snp_af_threshold \|\| SDP \< $depth_threshold\" --filterName \"SNP_AF\" -o $collect_snps_filtered\n";
		$string .= "gatk_toolkit -T VariantFiltration -R $reference -V $collect_indels --filterExpression \"SFREQ \< $indel_af_threshold \|\| SDP \< $depth_threshold\" --filterName \"INDEL_AF\" -o $collect_indels_filtered\n";
		# If finally it could be desirable to perform the final collecting without our inhouse scritp, comment next line (and uncomment last string line)
		$string .= "collecting_vcf_variants.py -i $collect_indels_filtered -i $collect_snps_filtered -o $final\n";
		#$string .= "gatk_toolkit -R $reference -T CombineVariants --variant:SNP $collect_snps_filtered --variant:INDEL $collect_indels_filtered -o $final -genotypeMergeOptions PRIORITIZE -priority INDEL,SNP\n";
	} else {
		$string .= "cp $collect $final\n";
	}
	my $job = launch_essay_job ($job_file, '', $previous_jobs, $string, $log_dir, $essay, $module);
	return [$job_file];
}

sub get_samparam {
	my ($essay,$readgroups,$bam) = @_;
	my $string;
	for (my $i = 0; $i < scalar(@{$readgroups}); $i++) {
		$string .= "<(samtools view -bh -L $essay->{target_reference}{target} -r $$readgroups[$i] $bam) ";
	}
	return $string;
}

1;