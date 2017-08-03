#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: intervals_coverage_depths.pl
# @Author: gmarco
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use json_handle;

my $config_json = $ENV{'MATRAZ_CONFIG'};

main();

sub main{
	
	#Test input args
	#-i 1.bam,2.bam,muestra_.bam -r test.bed -d 20 -threshold 65
	
	#Get input from menu. $beds is an array reference.
	my ($bams, $sample_name, $depth, $threshold, $panel, $user_csv, $output_dir, $num_samples) = menu();
	
	#Todo: Temporal fix, NC_stats module isn't ready for exomes atm.
	die 'ERROR: NC_stats module not ready for exome.' if ($panel =~ /exome/i);

	## MODULE: panel_info
	print "#############################################\n";
	print "START panel_info module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Loading panel, information: taget regions, CSV gene file, BED directories...\n";
	my ($target_bed, $csv, $gene_dir) = panel_info_json($panel);
	
	#If user defined CSV gene file list override config option:
	if ($user_csv ne ""){ 
		$csv = $user_csv;
	}
	
	print "DONE using panel_info: PANEL INFO LOADED OK.\n";
	print "END panel_info module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";

	print "#############################################\n";
	## Files & vars OK, Script starts here !
	print "INFO:\n";
	print "Using BAM file/s: @{$bams}\n";
	#print "Sample name: $sample_name\n";
	print "Panel: $panel\n";
	print "Depth value: $depth","x\n";
	print "Threshold: $threshold","\%\n";
	print "Output dir: $output_dir\n";
	
	## MODULE: samtools_depth
	print "#############################################\n";
	print "START samtools_depth module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Samtools Depth working...\n";
	
	my $output_samtools_depth = samtools_depth($bams, $target_bed, $sample_name, $depth, $threshold, $output_dir);
	
	print "DONE using samtools depth: $output_dir\/$output_samtools_depth\n";
	print "END samtools_depth module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	## MODULE: not_covered_bases
	print "#############################################\n";
	print "START not_covered_bases module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Calculating not covered bases...\n";
	
	my $output_not_covered_bases = not_covered_bases($output_samtools_depth, $sample_name, $depth, $threshold, $output_dir, $num_samples);
	
	print "DONE Calculating not covered bases: $output_dir\/$output_not_covered_bases\n";
	print "END covered_bases module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	## MODULE: intervals
	print "#############################################\n";
	print "START intervals module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Calculating intervals...\n";
	
	my $output_intervals = intervals($output_not_covered_bases, $sample_name, $depth, $threshold, $output_dir);
	
	print "DONE calculating intervals: $output_dir\/$output_intervals\n";
	print "END intervals module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	## MODULE: check_intervals
	check_intervals($output_intervals, $output_dir);
	
	## MODULE: output_hgmd	
	print "#############################################\n";
	print "START hgmd_main module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Obtaining HGMD information for non covered intervals...\n";

	my $output_hgmd = hgmd_main($output_intervals, $sample_name, $depth, $threshold, $output_dir);
	
	print "DONE getting HGMD information: $output_dir/$output_hgmd\n";
	print "END hgmd_main module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	## MODULE: chr_coverage
	print "#############################################\n";
	print "START chr_coverage: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Obtaining individual chr information for non covered intervals...\n";
	
	my $output_chr_coverage = chr_coverage($panel, $target_bed, $csv, $gene_dir, $output_intervals, $sample_name, $depth, $threshold, $output_dir);
	
	print "DONE getting chr information: $output_dir/$output_chr_coverage\n";
	print "END chr_coverage module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	#MODULE: vaporeta
	print "#############################################\n";
	print "START vaporeta module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Removing temporal files...\n";
	
	vaporeta($output_samtools_depth, $output_not_covered_bases, $output_dir);
	
	print "DONE removing temporal files.\n";
	print "END vaporeta module: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "#############################################\n";
	print "DONE: NON COVERED REGIONS SCRIPT\n";
}
	
sub menu{
	my (@bams, $sample_name, $depth, $threshold, $panel, $user_csv, $output_dir);
	
	$user_csv = "";
	
	#Get user input options
	GetOptions ("i=s" => \@bams, 'd=i' => \$depth,'threshold=i' => \$threshold, 'o=s' => \$output_dir, 
	'prefix=s' => \$sample_name, 'p=s' => \$panel, 'csv:s' => \$user_csv);
	#'r=s' => \$target_bed
	
	print "\n#############################################\n";
	print "######### SG Panel Non Covered Stats ########\n";
	print "#############################################\n";
	print "Remember that Depth & Threshold options minium value is 1.\n";
	print "NOTE: If just one sample is used, program will override user threshold and set it it 100%\n";
	print "If there's only one BAM file, program will check for each position if sample is covered or not with specified depth.\n";
	#print "\nFor each coverage score program checks it with the following expression: if (coverage >= depth) then it's well covered.\n";
	print "#############################################\n";
	print "\t",`date`,"\n";
	
	
	if (!@bams or !$sample_name or !$depth or !$threshold or !$output_dir or !$panel){
		print "-i BAM file/s, this input can be one or more files ie: -i file1.bam,file2.bam,file3.bam\n";
		print "-prefix Output_prefix ie: '-prefix prefix'\n";
		print "-p Panel name code: ie: '-p onco2'\n";
		print "-d Depth value to check coverage ie: for 10x '-d 10'\n";
		print "-threshold Coverage percentage threshold. ie: for 80% '-threshold 80'\n";
		print "-o output dir file ie: '-o mydir/output'\n\n";
		print "=============================================\n";
		print "Optional:\n";
		print "-csv csv gene file: '-csv my_gene_list.csv'\n\n";
		print "=============================================\n";
		print "Usage: sg_panel_not_covered_stats.pl -i 1.bam,2.bam,3bam -prefix sample_123 -p onco2 -d 20 -threshold 80 -o mydir/output\n\n";
		
		exit;
	}
	
	#Check BAM input
	if (!@bams){
		die "ERROR: You must specify at least one BAM input file !\n";
	}
	
	#Split multiple BAM list
	@bams = split(/,/,join(',',@bams));
	
	#Number of samples
	my $num_samples = scalar(@bams);
	
	#Check that BAM input files exists.
	foreach my $i (@bams){
		if (not -e $i){
			die "ERROR: Couldn't find BAM input file: $i\n";
		}
	}
	
	#Check output dir
	if (!defined $output_dir){
		die "ERROR: You must specify output directory !\n";
	}
	
	#Check output dir
	if (not -d $output_dir){
		print "WARNING: Couldn't find output directory !\n";
		print "Creating output directory: $output_dir\n";
		system ("mkdir -p $output_dir");
	
		if (not -d $output_dir){
			die "ERROR: Couldn't create output directory. Path incorrect or insufficient permissions.\n";
		}
	}
	
	if (!defined $depth){
		die "ERROR: You must specify the depth value !\n";
	}
	
	if (!defined $sample_name){
		die "ERROR: You must specify the sample name !\n";
	}
	
	if (!defined $threshold){
		die "ERROR: You must specify the threshold value !\n"
	}
	
	#Check that CSV exists
	if ($user_csv ne "" && not -e $user_csv){
		die "ERROR: Couldn't find CSV input file: $user_csv \n";
	}
	
	#Return values
	return (\@bams, $sample_name, $depth, $threshold, $panel, $user_csv, $output_dir, $num_samples);
	
}

sub vaporeta {
	my ($output_samtools_depth, $output_not_covered_bases, $output_dir) = @_;	  	
	#`rm $output_dir/$output_samtools_depth $output_dir/$output_not_covered_bases`;
	`rm $output_dir/$output_not_covered_bases`;
}

sub samtools_depth {
	my ($bams, $target_bed, $sample_name, $depth, $threshold, $output_dir) = @_;
	my $output_samtools_depth = "$sample_name"."_$depth"."x$threshold"."_samtools_depth.pileup";

	#`samtools depth -b $target_bed @{$bams} > $output_dir/$output_samtools_depth`;
	`samtools mpileup -A -d 100000 -l $target_bed @{$bams} > $output_dir/$output_samtools_depth`;
	return $output_samtools_depth;
}

sub not_covered_bases {
	#use feature qw(say);
	my ($output_samtools_depth, $sample_name, $depth, $threshold, $output_dir, $num_samples) = @_;
	
	my ($current_chr, $current_pos);
	my ($positive, $negative, $percentage);
	my ($col_number,$i);
	my @line;
	
	#my $output_not_covered_bases = "$sample_name"."_not_covered_bases_$depth"."x$threshold.pileup";
	my $output_not_covered_bases = "$sample_name"."_$depth"."x$threshold"."_not_covered_bases.pileup";
	
	#my $coverage = 1;
	
	#Open MPILUP input file to initialize vars.
	open MPILEUP, "$output_dir/$output_samtools_depth";
	open NOT_COVERED_BASES, ">$output_dir/$output_not_covered_bases";
	
	#Foreach line in MPILUP input
	while (<MPILEUP>) {
		chomp $_;
		
		#Store line info in @line array
		@line = split ("\t", $_);
		
		#Col number depends on the number of samples in input
		$col_number = scalar(@line);
		#$num_samples = $col_number-2;
		
		#If just one sample, threshold is 100%. So covered or not covered.
		if ($num_samples == 1){
			$threshold = 100;
			#print "Only one sample detected: THRESHOLD OVERWRITEN TO 100%\n";
		}
		
		#Current chr and position
		$current_chr = $line[0];
		$current_pos = $line[1];
		
		#Reset numeric vars
		($positive, $negative, $percentage)  = 0;
		
		#We check if the value is equal or higher to coverage for each sample on the line (position)
		for ($i = 3; $i < $col_number; $i=$i+3) {
			#print $line[$i],"<--------\n";
			if ($line[$i] >= $depth) {
				$positive++;
			}
 		}
 		
 		#We check if base must be kept or not
 		$percentage = ($positive/$num_samples)*100;
 		
 		#If percentage is below the threshold we report the non covered base position
 		#Use say because we don't want EOf end with \n empty line.
 		if ($percentage < $threshold){ 

 			print NOT_COVERED_BASES "$current_chr\t$current_pos\n";
 		}
	}
	
	#Close file handlers
	close MPILEUP;
	close NOT_COVERED_BASES;
	
	return $output_not_covered_bases;

}

sub intervals {
	my ($output_not_covered_bases, $sample_name, $depth, $threshold, $output_dir) = @_;
	
	my ($current_chr, $start_chr, $prev_chr);
	my ($current_pos, $start_pos, $prev_pos);
	my @line;

	my $output_intervals = "$sample_name"."_$depth"."x$threshold"."_not_covered_intervals.bed";
	
	
	open NOT_COVERED_BASES, "$output_dir/$output_not_covered_bases";
	open NOT_COVERED_INTERVALS, ">$output_dir/$output_intervals";
	while(<NOT_COVERED_BASES>){
		chomp $_;
		
		$start_chr = $current_chr = $prev_chr = (split ("\t",$_))[0];
		$start_pos = $current_pos = $prev_pos = (split ("\t",$_))[1];
		
		chomp ($start_chr, $start_pos);
		
		last;
	}
	
	#We come back to file start
	seek NOT_COVERED_BASES, 0, 0;
	
	#Foreach line in MPILUP input
	while (<NOT_COVERED_BASES>) {
		chomp $_;
		
		@line = split ("\t", $_);
		
		#Update prev vars
		$prev_pos = $current_pos;
		$prev_chr = $current_chr;
		
		#Update current vars
		$current_chr = $line[0];
		$current_pos = $line[1];
		
		chomp ($current_chr, $current_pos);
		
		#If there's a gap greater than 1 between positions we create interval
		if (($current_pos - $prev_pos) > 1){
			print NOT_COVERED_INTERVALS "$start_chr\t$start_pos\t$prev_pos\n";
			$prev_chr = $current_chr;
			$start_pos = $current_pos;
			$start_chr = $current_chr;
		}
		 
		
		#If there's a chromosome change we create interval
		if ($prev_chr ne $current_chr){
			print NOT_COVERED_INTERVALS "$start_chr\t$start_pos\t$prev_pos\n";
			$prev_chr = $current_chr;
			$start_pos = $current_pos;
			$start_chr = $current_chr;
		}
		
		#If EOF create remaining interval
		if (eof(NOT_COVERED_BASES)){
			print NOT_COVERED_INTERVALS "$start_chr\t$start_pos\t$current_pos\n";
		} 
	}
	
	#Close file handlers
	close NOT_COVERED_BASES;
	close NOT_COVERED_INTERVALS;
	
	return $output_intervals;
}

sub check_intervals {
	my ($output_intervals, $output_dir) = @_;
	
	if (-z "$output_dir/$output_intervals"){
		print "#############################################\n";
		print "$output_intervals file is empty !\n";
		print "All bases have enough depth. You may be using a very low depth value !\n";
		exit;
	}
}
	
sub hgmd_main {
	
	
	
	my ($output_intervals, $sample_name, $depth, $threshold, $output_dir) = @_;
	#my $masterbiobed = "/share/apps/scripts/genomics/panels/intervals_coverage_depths/masterbioBed_sorted.bed";
	my $masterbiobed = $ENV{'MASTERBIOBED'};
	my $hgmd_not_covered_mutations = "$sample_name"."_$depth"."x$threshold"."_not_covered_hgmd_mutations.bed";
	my $intervals_not_covered_bases = "$sample_name"."_$depth"."x$threshold"."_not_covered_intervals_bases_column.bed";
	my $intervals_not_covered_merged_bases = "$sample_name"."_$depth"."x$threshold"."_not_covered_merged_bases.bed";
	my $output_hgmd = "$sample_name"."_$depth"."x$threshold"."_not_covered_hgmd_final.bed";

	#Obtain HGMD information for non covered regions...
	`intersectBed -a $output_dir/$output_intervals -b $masterbiobed -wa -wb | awk 'BEGIN { OFS = "\\t";} {print \$1,\$2,\$3,\$7,\$8,\$9,\$10,\$11}' > $output_dir/$hgmd_not_covered_mutations`;
	
	print "DONE ! HGMD information for non covered design interval: $output_dir/$hgmd_not_covered_mutations\n";
	print "=============================================\n";
	print "Obtaining number of bases for each non covered interval...\n";
	
	#Obtain number of bases in column format for non_covered regions.
	system("awk '{print \$3-\$2+1}' $output_dir/$output_intervals > $output_dir/$intervals_not_covered_bases");
	
	print "DONE ! Number of bases for each non covered interval: $output_dir/$intervals_not_covered_bases\n";
	print "=============================================\n";
	
	#Merge non_covered regions with bases count for each line
	`paste -d'\\t' $output_dir/$output_intervals $output_dir/$intervals_not_covered_bases > $output_dir/$intervals_not_covered_merged_bases`;
	
	print "Merging final HGMD resulsts file...\n";
	hgmd_final_merge($intervals_not_covered_merged_bases, $hgmd_not_covered_mutations, $output_hgmd, $output_dir);
	
	#Remove tmp hgmd_main files 
	`rm $output_dir/$intervals_not_covered_bases $output_dir/$intervals_not_covered_merged_bases $output_dir/$hgmd_not_covered_mutations`;
	
	print "DONE hgmd_final_merge !\n";
	print "=============================================\n";
	
	return $output_hgmd;
}

sub hgmd_final_merge {
	my ($intervals_not_covered_merged_bases, $hgmd_not_covered_mutations, $output_hgmd, $output_dir) = @_;
	
	open BASE, "$output_dir/$intervals_not_covered_merged_bases";
	open HGMD_FINAL_OUTPUT, ">$output_dir/$output_hgmd";
	print HGMD_FINAL_OUTPUT "#CHR\tSTART_POS\tEND_POS\tNUM_BASES\tHGMD\n";
	
	while (<BASE>){
		chomp $_;
		my $found = 0;
		my @split_base = split("\t", $_);

		open HGMD, "$output_dir/$hgmd_not_covered_mutations";
		
		while (<HGMD>){
			chomp $_;
			my @split_hgmd = split("\t", $_);
			if ($split_base[0] eq $split_hgmd[0] && $split_base[1] == $split_hgmd[1] && $split_base[2] == $split_hgmd[2]){
				$found = 1;
				
				my ($count, $hgmd_desc_size);
				my $hgmd_desc = "";
				$hgmd_desc_size = scalar (@split_hgmd);
				
				#Iterate over all HGMD description
				if ($split_hgmd[3]){
					for ($count = 3; $count < $hgmd_desc_size; $count++) {
	 					$hgmd_desc .= " $split_hgmd[$count]";
 					}
				}
 				
				#print HGMD_FINAL_OUTPUT "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\t$split_hgmd[3]\n";
				if ($hgmd_desc ne ""){
					print HGMD_FINAL_OUTPUT "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\t$hgmd_desc\n";
				}
				else {
					print HGMD_FINAL_OUTPUT "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\n";
				}
			} 	
		}
		close HGMD;
		
		if ($found == 0) {
			print HGMD_FINAL_OUTPUT "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\n";
		}	
	}
	close HGMD;
	close BASE;
	close HGMD_FINAL_OUTPUT;
}

#Todo: Depecrated function, panel_info_json reads config from Matraz's config JSON.
#sub panel_info {
#	my $panel = lc($_[0]);
#	my (@line, $target_bed, $csv, $gene_dir);
#	my @panels;
#	
#	open PANEL_DB, "/share/apps/scripts/genomics/panels/intervals_coverage_depths/panel_db.txt";
#	#open PANEL_DB, "panel_db.txt";
#	
#	while (<PANEL_DB>){
#		if ($_ =~ /^#/) { next; }
#		else {
#			 @line = split (/\|/,$_);
#			 
#			 #Store available panel names to print if user is wrong.
#			 push(@panels, $line[0]);
#
#			 if ($line[0] eq $panel){ 
#			 	$target_bed = $line[1];
#			 	$csv = $line[2];
#			 	$gene_dir = $line[3];
#			 	
#			 	close PANEL_DB;
#			 	return ($target_bed, $csv, $gene_dir);
#			 }
#		}
#	}
#	die "ERROR: Panel $panel doesn't exist or is not ready to be used yet !\nCurrent available panels: @panels\n\n";
#}

sub panel_info_json {
	my $panel = lc($_[0]);
	
	my %analysis_type = get_from_json($config_json);
	my (@available_list, $available_list);
	
	#Obtain a list of available target regions in case user is not selecting any valid option.
	foreach my $key (keys %analysis_type ) {
		foreach my $target_available (sort keys %{$analysis_type{$key}}) {
			push (@available_list, $target_available);
		}
	}
	
	$available_list = join(", ",@available_list);
	
	die "ERROR: $panel doesn't exist in JSON config file !\nAvailable analisis list: $available_list.\n" if not exists ${analysis_type{'target_reference'}{$panel}};
	die "ERROR: $panel CSV file not defined or empty in Matraz config !\n" if not exists ${analysis_type{'target_reference'}{$panel}{genes_csv}};
	die "ERROR: $panel BED gene dir not defined or empty in Matraz config !\n" if not exists ${analysis_type{'target_reference'}{$panel}{bed_gene_dir}};
	
	#Get required information panel information from Matraz config JSON
	my $target_bed = ${analysis_type{'target_reference'}{$panel}{target}};
	my $csv = ${analysis_type{'target_reference'}{$panel}{genes_csv}};
	my $gene_dir = ${analysis_type{'target_reference'}{$panel}{bed_gene_dir}};
	
	return ($target_bed, $csv, $gene_dir);
	
}

sub chr_coverage {
	
	#Input vars
	my ($panel, $target_bed, $csv, $gene_dir, $output_intervals, $sample_name, $depth, $threshold, $output_dir) = @_;
	
	#Vars
	my $gene;
	my ($chr,$input_totals_coords,$input_exon_coords,$input_splicing_coords,$input_strong_splicing_coords,
	$input_UTR_coords,$input_vars_coords, $input_protein_coding_exon);
	
	my $output_chr_coverage = "$sample_name"."_$depth"."x$threshold"."_chr_not_covered_stats.txt";
	my $output_exon_coverage = "$sample_name"."_$depth"."x$threshold"."_exons_not_covered_stats.txt";
	
	#Gene list for panel
	my @gene_list = &get_gene_list($csv);
	
	chomp $gene_dir;
	
	check_existing_gene_dirs($gene_dir, @gene_list);

	my $awk_command = " | awk '{print \$3-\$2+1}'";
	my $awk = "awk '{print \$3-\$2+1}'";
	
	## Open output_file
	open CHR_COVERAGE_OUTPUT, ">$output_dir/$output_chr_coverage";
	open EXON_COVERAGE_OUTPUT, ">$output_dir/$output_exon_coverage";
	print CHR_COVERAGE_OUTPUT "#HGNC\tCHR\tTotal_Non_Covered_Gene\tTotal_Gene\tPercentage_Non_Covered_Gene\tTotal_Non_Covered_Exon\tTotal_Exon\tPercentage_Non_Covered_Exon\tTotal_Non_Covered_Protein_Coding_Exon\tTotal_Protein_Cooding_Exon\tPercentage_Non_Covered_Protein_Coding_Exon\tTotal_Non_Covered_Splicing\tTotal_Splicing\tPercentage_Non_Covered_Splicing\tTotal_Non_Covered_Strong_Splicing\tTotal_Strong_Splicing\tPercentage_Non_Covered_Strong_Splicing\tTotal_Non_Covered_UTR\tTotal_UTR\tPercentage_Non_Covered_UTR\n";
	print EXON_COVERAGE_OUTPUT "#HGNC\tCHR\tTranscript_name\tExon_name\tTotal_Non_Covered_Exon\tTotal_Exon\tPercentage_Non_Covered_Exon\n"; 
	
	
	foreach $gene (@gene_list){
		
		chomp $gene;		
		
		$input_totals_coords = $gene_dir."/$gene/target.bed";
		$input_exon_coords = $gene_dir."/$gene/exon.bed";
		$input_protein_coding_exon = $gene_dir."/$gene/cds.bed"; 
		$input_splicing_coords = $gene_dir."/$gene/splicing.bed";
		$input_strong_splicing_coords = $gene_dir."/$gene/strong_splicing.bed";
		$input_UTR_coords = $gene_dir."/$gene/UTR.bed";

		$chr = (split ("\n",`awk '{print \$1}' $input_totals_coords`))[0];
 		
 		#Get exon file bed list
 		my $transcript_dir = $gene_dir."/$gene/canonical_transcript";
 		my @transcript_list = get_file_dir_list($transcript_dir);
 		#my @exon_files = get_file_dir_list($exon_dir);
 		
 		#For each exon we calculate all stuff
 		foreach my $transcript (@transcript_list){
 			print "Processing transcript: $transcript for gene: $gene\n";
 			my $exon_dir = $gene_dir."/$gene/canonical_transcript/$transcript/exons";
 			my @exon_files = get_file_dir_list($exon_dir);
 			
 			foreach my $single_exon_file (@exon_files){
 				my $full_path_exon_file = $exon_dir."/$single_exon_file";
 				
 				if (-z $full_path_exon_file){
 					next;
 				}
 				
 				(my $exon_basename, my $dir, my $ext) = fileparse($single_exon_file, qr/\.[^.]*/);
 				my $sum_single_exon = col_sum(`$awk $full_path_exon_file`);
 				my $sub_single_exon = col_sum(`subtractBed -a $full_path_exon_file -b $output_dir/$output_intervals $awk_command`);
				#print "subtractBed -a $full_path_exon_file -b $output_dir/$output_intervals $awk_command";
 				my ($perc_single_exon, $resta_single_exon) = &percentage_coverage($sum_single_exon, $sub_single_exon);
 				
 				#If number of non covered exon bases is different of 0 then report it.
 				if ($resta_single_exon != 0){
 					print EXON_COVERAGE_OUTPUT "$gene\t$chr\t$transcript\t$exon_basename\t$resta_single_exon\t$sum_single_exon\t$perc_single_exon\n"; 					
 				}
 			}
 		}
 		
 		#Chr non covered stats
		my ($sub_total_gen, $sub_total_exon_coords, $sub_total_protein_coding_exon_coords, $sub_total_splicing_coords, $sub_total_strong_splicing_coords, $sub_total_UTR_coords) = (0, 0 ,0 ,0 ,0 ,0);
		my $sum_total_gen = col_sum(`$awk $input_totals_coords`);
		my $sum_total_exon_coords = col_sum(`$awk $input_exon_coords`);
		my $sum_total_protein_coding_exon_coords = col_sum(`$awk $input_protein_coding_exon`);
		my $sum_total_splicing_coords = col_sum(`$awk $input_splicing_coords`);
		my $sum_total_strong_splicing_coords = col_sum(`$awk $input_strong_splicing_coords`);
		my $sum_total_UTR_coords = col_sum(`$awk $input_UTR_coords`);
		
		
		#Check for each input intersect file that BED file is not empty.
		if (!-z $input_totals_coords){
			$sub_total_gen = col_sum(`subtractBed -a $input_totals_coords -b $output_dir/$output_intervals $awk_command`);	
		}
		
		if (!-z $input_exon_coords){
			$sub_total_exon_coords = col_sum(`subtractBed -a $input_exon_coords -b $output_dir/$output_intervals $awk_command`);
		}
		
		if (!-z $input_protein_coding_exon){
			$sub_total_protein_coding_exon_coords = col_sum(`subtractBed -a $input_protein_coding_exon -b $output_dir/$output_intervals $awk_command`);
		}
		
		if (!-z $input_splicing_coords){
			$sub_total_splicing_coords = col_sum(`subtractBed -a $input_splicing_coords -b $output_dir/$output_intervals $awk_command`);
		}
		
		if (!-z $input_strong_splicing_coords){
			$sub_total_strong_splicing_coords = col_sum(`subtractBed -a $input_strong_splicing_coords -b $output_dir/$output_intervals $awk_command`);
		}
		
		if (!-z $input_UTR_coords){
			$sub_total_UTR_coords = col_sum(`subtractBed -a $input_UTR_coords -b $output_dir/$output_intervals $awk_command`);
		} 

		my ($perc_gen, $resta_gen) = &percentage_coverage($sum_total_gen, $sub_total_gen);
		my ($perc_exon, $resta_exon) = &percentage_coverage($sum_total_exon_coords, $sub_total_exon_coords);
		my ($perc_protein_coding_exon, $resta_protein_coding_exon) = &percentage_coverage($sum_total_protein_coding_exon_coords, $sub_total_protein_coding_exon_coords);
		my ($perc_splicing, $resta_splicing) = &percentage_coverage($sum_total_splicing_coords, $sub_total_splicing_coords);
		my ($perc_strong_splicing, $resta_strong_splicing) = &percentage_coverage($sum_total_strong_splicing_coords, $sub_total_strong_splicing_coords);
		my ($perc_UTR, $resta_UTR) = &percentage_coverage($sum_total_UTR_coords, $sub_total_UTR_coords);
		
		print CHR_COVERAGE_OUTPUT "$gene\t$chr\t$resta_gen\t$sum_total_gen\t$perc_gen\t$resta_exon\t$sum_total_exon_coords\t$perc_exon\t$resta_protein_coding_exon\t$sum_total_protein_coding_exon_coords\t$perc_protein_coding_exon\t$resta_splicing\t$sum_total_splicing_coords\t$perc_splicing\t$resta_strong_splicing\t$sum_total_strong_splicing_coords\t$perc_strong_splicing\t$resta_UTR\t$sum_total_UTR_coords\t$perc_UTR\n";
	}
	#Close file handlers
	close CHR_COVERAGE_OUTPUT;
	close EXON_COVERAGE_OUTPUT;
	return $output_chr_coverage;
}

sub get_file_dir_list {
	my $dir = shift;
	my @dir_file_list;
	
	opendir (DIR, $dir) or die "ERROR: $dir not found, you're probably missing gene bed directories..\n";
	
	while (my $file = readdir(DIR)) {
		
		#Skip . and .. special folders
		if ($file eq "." or $file eq ".."){
	 		next;
	 	}
	 	
	 	push (@dir_file_list, $file);
	 	
	}
	return @dir_file_list;
}

sub check_existing_gene_dirs {
	my ($gene_dir, @gene_list) = @_;
	
	foreach my $gene (@gene_list){
		my $gene_dir_path = $gene_dir."/$gene";
		is_folder_empty($gene_dir_path);
	}
}

sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "ERROR: $dirname is not a directory or not exsists !!\n";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

sub percentage_coverage {
	my ($total_elemento, $total_subtract) = @_;
	my $resta = $total_elemento - $total_subtract;
	
	if ($resta == 0){
		return 0.00, $resta
	}
	
	my $porc = ($resta/$total_elemento)*100;
	$porc = sprintf("%.2f", $porc);
	return $porc, $resta;
}

sub get_gene_list {
	my $csv = $_[0];
	my @gene_list;
	
	open GENE_CSV, "$csv" or die "ERROR: Couldn't open CSV gene file !\n";
	while (<GENE_CSV>){
		chomp $_;
		my @gene_split = split(',', $_);
		push (@gene_list, $gene_split[0]); 		
	}
	close GENE_CSV;
	
		
	return @gene_list;
}

sub col_sum {
	my @columna = @_;
	my $valor = 0;
	my $sumatorio = 0;
	my $elemento; 
	
	foreach $elemento (@columna){
		chomp $elemento;
		$valor = $elemento;
		$sumatorio = $sumatorio + $valor;	
	}
	return $sumatorio;
}
