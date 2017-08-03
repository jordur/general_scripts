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
my ($rawmetrics,$merged,$pcr);
my $output = cwd() . "/output_file.tsv";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"r=s"	=> \$rawmetrics,
				"m=s"	=> \$merged,
				"pcr=s"	=> \$pcr,
				"o=s"	=> \$output,
);

print "\n=================================================================
$basename: Script for obtaining sequencing lane metrics v0.1\n";

if (not defined($rawmetrics) or not defined($merged) or not defined($pcr)) {
	die "
			Parameters:
		-r    rawdata metrics files (comma separated list of outputs of lane_checker.pl)
		-m    merged bam file of all lanes belonging to the replicate
		-pcr  merged bam file of all lanes belonging to the replicate after performing removal of PCR duplicates
		-o    output name (optional)
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
	print $out "Total_reads\tMapped_reads\tPercentage_of_mapped_reads\tReads_after_duplicate_removal\tPercentage_of_reads_after_duplicate_removal\n";
	
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
	my $pcr_reads = `bamtools stats -in $pcr | grep "^Mapped" | awk '{print \$3}'`;
	chomp($pcr_reads);
	
	# Print mapping metrics to output
	my $percentage_merged = sprintf "%.2f",$merged_reads/$replicate_reads*100;
	my $percentage_pcr = sprintf "%.2f",$pcr_reads/$replicate_reads*100;
	print $out $replicate_reads."\t".$merged_reads."\t".$percentage_merged."\t".$pcr_reads."\t".$percentage_pcr."\n";
		
	# Close output file handle
	$out->close();
}

sub GetMetrics(){
	my ($file,$type) = @_;
	my %metrics;
	
	# Open file for reading
	my $fh = &get_in_file_handle($file);
	
	# Check each line
	my $reads = 0; # Number of reads in file
	my $nts = 0; # Number of nucleotides per read
	my $current_read;
	while (<$fh>){
		my $line = $_;
		chomp ($line);
		if (substr($line,0,1) eq ">"){
			$current_read = $line;
			$reads++;
		} else {
			my $current_nts;
			if ($type eq "fasta"){
				$current_nts = length($line)-1;
			} elsif ($type eq "qual"){
				my @quals = split(" ",$line);
				$current_nts = $#quals+1;
			}
			if ($nts == 0){
				$nts = $current_nts;
			} elsif ($current_nts != $nts){
				die "ERROR: Invalid number of nucleotides in file $file in read $current_read (".length($line).")!!\n";
			}
		}
	}
	$fh->close();
	
	$metrics{read_length} = $nts;
	$metrics{reads} = $reads;
	return %metrics;
}

exit;
