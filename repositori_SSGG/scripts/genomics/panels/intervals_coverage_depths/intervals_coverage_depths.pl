#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: intervals_coverage_depths.pl
# @Author: Guille
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use Getopt::Long;
use POSIX qw(strftime);

main();
#menu();
#covered_bases();
#intervals();


sub main{
	
	#Test input args
	#-i 1.bam,2.bam,muestra_.bam -r test.bed -d 20 -threshold 65
	
	#Get input from menu. $beds is an array reference.
	my ($bams, $target_bed, $sample_name, $depth, $threshold, $output_dir) = menu();
	
	print "Using BAM file/s: @{$bams}\n";
	print "Target region BED: $target_bed\n";
	print "Sample name: $sample_name\n";
	print "Depth value: $depth","x\n";
	print "Threshold: $threshold","\%\n";
	print "Output dir: $output_dir\n";
	
	#Call samtools depth to get the bed with depth values.	
	my $output_samtools_depth = samtools_depth($bams, $target_bed, $sample_name, $output_dir);
	my $output_not_covered_bases = not_covered_bases($output_samtools_depth, $sample_name, $depth, $threshold, $output_dir);
	my $output_intervals = intervals($output_not_covered_bases, $sample_name, $depth, $threshold, $output_dir);
	

sub menu{
	my (@bams, $target_bed, $sample_name, $depth, $threshold, $output_dir);
	
	#Get user input options
	GetOptions ("i=s" => \@bams, 'r=s' => \$target_bed,
	'd=i' => \$depth,'threshold=i' => \$threshold, 'o=s' => \$output_dir, 
	's=s' => \$sample_name);
	
	
	print "#############################################\n";
	print "######### Intervals Depth Coverage ##########\n";
	print "#############################################\n";
	print "\t",`date`,"\n";
	
	if (!@bams or !$target_bed or !$sample_name or !$depth or !$threshold or !$output_dir){
		print "-i BAM file/s, this input can be one or more files ie: -i file1.bam,file2.bam,file3.bam\n";
		print "-r BED target regions file\n";
		print "-s Sample name ie: Oto2_1234\n";
		print "-d Depth value to check coverage ie: for 10x '-d 10'\n";
		print "-threshold Coverage percentage threshold. ie: for 80% '-c 80'\n";
		print "-o output dir file ie: /share/gluster/OnGoing/BF123_example/coverage_stats\n";
		print "#############################################\n";
		exit;
	}
	
	#Check BAM input
	if (!@bams){
		die "ERROR: You must specify at least one BAM input file !\n";
	}
	
	#Check output dir
	if (!defined $output_dir){
		die "ERROR; You must specify output directory !\n";
	}
	
	#Check output dir
	if (not -d $output_dir){
		print "WARNING: Couldn't find output directory !\n";
		print "Creating output directory: $output_dir\n";
		print "#############################################\n";
		system ("mkdir -p $output_dir");
	
		if (not -d $output_dir){
			die "ERROR: Couldn't create output directory. Path incorrect or insufficient permissions.\n";
		}
	}
	
	#Check all input values	
	if (!defined $target_bed){
		die "ERROR: You must specify BED target region file !\n";
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
	
	#Split multiple BAM list
	@bams = split(/,/,join(',',@bams));
	
	#Check that BAM input files exists.
	foreach my $i (@bams){
		if (not -e $i){
			die "ERROR: Couldn't find BAM input file: $i\n";
		}
	}
	
	#Check that BED target regions file exists.
	if (not -e $target_bed){
		die "ERROR: Couldn't find BED target regions file: $target_bed\n";
	}
	
	#Return values
	return (\@bams, $target_bed, $sample_name, $depth, $threshold, $output_dir);
	
}

sub samtools_depth {
	my ($bams, $target_bed, $sample_name, $output_dir) = @_;
	my $output_samtools_depth = "$sample_name"."_samtools_depth.pileup";
	
	print "#############################################\n";
	print "START samtools_depth: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Samtools Depth working...\n";
	`samtools depth -b $target_bed @{$bams} > $output_dir/$output_samtools_depth`;
	print "DONE Samtools Depth: $output_dir\/$output_samtools_depth\n";
	print "END samtools_depth: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	return $output_samtools_depth;
}

sub not_covered_bases {
	#use feature qw(say);
	my ($output_samtools_depth, $sample_name, $depth, $threshold, $output_dir) = @_;
	
	my ($current_chr, $current_pos);
	my ($positive, $negative, $percentage);
	my ($col_number,$num_samples,$i);
	my @line;
	
	#my $output_not_covered_bases = "$sample_name"."_not_covered_bases_$depth"."x$threshold.pileup";
	my $output_not_covered_bases = "$sample_name"."_$depth"."x$threshold"."_not_covered_bases.pileup";
	
	#my $coverage = 1;

	print "#############################################\n";
	print "START not_covered_bases: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Calculating not covered bases for sample: $sample_name\n";
	
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
		$num_samples = $col_number-2;
		
		#Current chr and position
		$current_chr = $line[0];
		$current_pos = $line[1];
		
		#Reset numeric vars
		($positive, $negative, $percentage)  = 0;
		
		#We check if the value is equal or higher to coverage for each sample on the line (position)
		for ($i = 2; $i < $col_number; $i++) {
			if ($line[$i] >= $depth) {
				$positive++;
			}
			
 		}
 		
 		#We check if base must be kept or not
 		$percentage = ($positive/$num_samples)*100;
 		
 		#If percentage is equal below the threshold we report the non covered base position
 		#Use say because we don't want EOf end with \n empty line.
 		if ($percentage < $threshold){ 
 			print NOT_COVERED_BASES "$current_chr\t$current_pos\n";
 		}
	}
	
	print "DONE Calculating not covered bases: $output_dir\/$output_not_covered_bases\n";
	print "END covered_bases: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
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
	#my ($col_number,$num_samples,$i);

	#my ($positive, $negative, $percentage);
	#my $coverage = 1;

	my $output_intervals = "$sample_name"."_$depth"."x$threshold"."_not_covered_intervals.bed";
	
	print "#############################################\n";
	print "START intervals: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	print "Calculating intervals...\n";	
	
	#
	open NOT_COVERED_BASES, "$output_dir/$output_not_covered_bases";
	open NOT_COVERED_INTERVALS, ">$output_dir/$output_intervals";
	while(<NOT_COVERED_BASES>){
		chomp $_;
		
		#Col number depends on the number of samples in input
		#$col_number = scalar(split ("\t", $_));
		#$num_samples = $col_number-2;
		
		$start_chr = $current_chr = $prev_chr = (split ("\t",$_))[0];
		$start_pos = $current_pos = $prev_pos = (split ("\t",$_))[1];
		
		last;
		#exit;
	}
	
	#We come back to file start
	seek NOT_COVERED_BASES, 0, 0;
	
	#Foreach line in MPILUP input
	while (<NOT_COVERED_BASES>) {
		chomp $_;
		
		#Reset numeric vars
		#($positive, $negative, $percentage)  = 0;
		
		@line = split ("\t", $_);
#		
#		#We check if the value is equal or higher to coverage for each sample on the line (position)
#		for ($i = 2; $i < $col_number; $i++) {
#			#print "$line[$i] ";
#			if ($line[$i] >= $coverage) {
#				$positive++;
#			}
# 		}
# 		
# 		$percentage = ($positive/$num_samples)*100;
# 		if ($percentage >= 50){
# 			#print "\t$positive\/$num_samples\t$percentage%\n";
# 		}
		
		
		#Update prev vars
		$prev_pos = $current_pos;
		$prev_chr = $current_chr;
		
		#Update current vars
		$current_chr = $line[0];
		$current_pos = $line[1];
		#$current_chr = (split ("\t",$_))[0];
		#$current_pos = (split ("\t",$_))[1];
		
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
	
	print "DONE intervals: $output_dir\/$output_intervals\n";
	print "END intervals: ",strftime("%a, %d\/%m\/%Y %H:%M:%S %z", localtime(time()))."\n";
	
	return $output_intervals;
	}
}