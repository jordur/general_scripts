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
use Math::Round;
use file_handle;

# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my ($rawmetrics,$merged,$pcr,$q1,$target);
my $output = cwd() . "/output_file.tsv";
my $tmp_folder = cwd();

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"r=s"	=> \$rawmetrics,
				"m=s"	=> \$merged,
				"q1=s"	=> \$q1,
				"pcr=s"	=> \$pcr,
				"t=s"	=> \$target,
				"o=s"	=> \$output,
				"f=s"	=> \$tmp_folder,
);

print "\n=================================================================
$basename: Script for obtaining sequencing lane metrics v0.1\n";

if (not defined($rawmetrics) or not defined($merged) or not defined($pcr) or not defined($target)) {
	die "
			Parameters:
		-r    rawdata metrics files (comma separated list of outputs of lane_checker.pl)
		-m    merged bam file of all lanes belonging to the replicate
		-q1   merged bam file of all lanes belonging to the replicate after performing removal of reads with low qualities (Q1). Optional parameter
		-pcr  merged bam file of all lanes belonging to the replicate after performing removal of PCR duplicates (and realignment+recalibration)
		-t    target regions bed file
		-o    output name (optional)
		-f    temp folder (optional)
\n".localtime()."\n=================================================================\n\n";
}

&main ();
print  "\nINFO: Process successfully finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	my $path = `pwd`;
	chomp $path;
	
	# Print output header
	my $out = &get_out_file_handle($output);
	print $out "Total reads\tMapped reads\tPercentage of mapped reads\tReads after low quality reads removal\tPercentage of reads after low quality reads removal\tReads after duplicate removal\tPercentage of reads after duplicate removal\tReads on target\tPercentage of reads on target\n";
	
	# Get rawdata metrics
	my @rawdata_metrics = split(",",$rawmetrics);
	my $replicate_reads = 0;
	foreach my $file (@rawdata_metrics){
		my $fh = get_in_file_handle ($file);
		
		my $line = 0;
		my $lane_reads = 0;
		while (<$fh>){
			my $file_line = $_;
			chomp($file_line);
			$line ++;
			# Second and fourth line of file must contain the number of F3 and F5 reads
			if ($line == 2 || $line == 4){
				my @tmp = split('\t',$file_line);
				$lane_reads = $lane_reads + $tmp[2];
			}
		}
		$fh->close();
		
		$replicate_reads = $replicate_reads + $lane_reads;
	}
	
	# Get bam metrics
	my $merged_reads = `bamtools stats -in $merged | grep "^Mapped" | awk '{print \$3}'`;
	chomp($merged_reads);
	my $q1_reads = "NULL";
	if (defined($q1)){
		$q1_reads = `bamtools stats -in $q1 | grep "^Mapped" | awk '{print \$3}'`;
		chomp($q1_reads);
	}
	my $pcr_reads = `bamtools stats -in $pcr | grep "^Mapped" | awk '{print \$3}'`;
	chomp($pcr_reads);
	
	#my $log = `sg_reads_on_target_bamtools.pl -p $target -i $pcr -o $tmp_folder/on_target -t 1`;
	my $log = `reads_on_target_bedtools.py -p $target -i $pcr -o $tmp_folder/on_target.stat.all -f $tmp_folder`;
	my $on_target_reads = `awk '{suma+=\$2} END {print suma}' $tmp_folder/on_target.stat.all`;
	chomp($on_target_reads);
	
	# Print mapping metrics to output
	my $percentage_merged = sprintf "%.2f",$merged_reads/$replicate_reads*100;
	my $percentage_q1 = "NULL";
	if (defined($q1)){
		$percentage_q1 = sprintf "%.2f",$q1_reads/$replicate_reads*100;
	}
	my $percentage_pcr = sprintf "%.2f",$pcr_reads/$replicate_reads*100;
	my $on_target = sprintf "%.2f",$on_target_reads/$replicate_reads*100;
	print $out $replicate_reads."\t".$merged_reads."\t".$percentage_merged."\t".$q1_reads."\t".$percentage_q1."\t".$pcr_reads."\t".$percentage_pcr."\t".$on_target_reads."\t".$on_target."\n";
		
	# Close output file handle
	$out->close();
}

exit;