#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script for getting and collecting BioScope indels from different samples in target regions
# @Author: Arbol
# @Contributors: JM Rosa
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
use utils::bioinfo_files_utilities qw(import_intervals_file create_intervals_file);
use List::Util qw[min max];

my $intervals_file;
my $vars_file;
my $output_id;

GetOptions (
				"f=s"	=> \$intervals_file,
				"v=s"	=> \$vars_file,
				"o=s"	=> \$output_id,
);

print "\n\n=================================================================
Getting and collecting BioScope indels from different samples in target regions v1.0\n";

if (not defined $intervals_file or not defined $vars_file) {
	die "
			Options:
		-f chr-coord file (chrXX:coord1-coord2 files)
		-v vars file (file containing the called variants)
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;

if (not defined $output_id) {
	$output_id = "found_variants.tsv";
}

my %control_intervals;
my $log_file = "Log.txt";
open my $LOG, '>', $log_file or die ("ERROR: Impossible to create $log_file!! $!");
open my $OUT, '>', $output_id or die ("ERROR: Impossible to create $vars_file!! $!");

import_intervals_file($intervals_file, "chr-coord", \%control_intervals, $LOG);

my @chroms = (1..22, "X", "Y", "M");
foreach ( @chroms ) {
	my $current_chrom = $_;
	if ( exists($control_intervals{$current_chrom}) ) {
		for (my $i=0;$i<=$#{$control_intervals{$current_chrom}};$i++){
			my $begin = min(${$control_intervals{$current_chrom}}[$i][0],${$control_intervals{$current_chrom}}[$i][1]);
			my $end = max(${$control_intervals{$current_chrom}}[$i][0],${$control_intervals{$current_chrom}}[$i][1]);
			print "chr$current_chrom : $begin - $end \n";
			for (my $j=$begin;$j<=$end;$j++){
				#print "grep chr$current_chrom $vars_file | grep $j\n";
				my $vars = `grep chr$current_chrom $vars_file | grep $j`;
				print $OUT $vars;
			}
		}
	}
}

exit;

# Open log, output and error files
print LOG "sg_parsing_gatk.pl -v -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

# First, get all sample names and set file pointers to each file
my @samples; # Array containing the sample names of the files
my @files; # Array containing the handles to the vcf files
my @lines; # Array containing one line of each vcf file
my @variant; # Array containing the corresponding variant
my $SB; #Strand bias value for each indel
my $samples_for_indel; #Number of samples with a given indel


# Print output header
print OUT "#Chr\tPosition\tRef_Allele\tVar_Allele\tType";
foreach my $sample (@samples) {
	print OUT "\t$sample\_depth\t$sample\_geno\t$sample\_freq";
}
print OUT "\tSB\n";

# Now parse all sample files and output indels
my $first_pos;
my $first_variant;
do {
	$first_variant = variants->get_first_variant(\@variant);
	if (defined($first_variant)){
		my $line = "chr$first_variant->{chr}\t$first_variant->{position}\t$first_variant->{ref_allele}\t$first_variant->{var_allele}\t-";
		$SB = 0;
		$samples_for_indel = 0;
		for (my $i=0;$i<=$#lines;$i++){
			if (defined($variant[$i])){
				if (variants->compare($variant[$i],$first_variant) == 0){
					$line = $line . "\t$variant[$i]->{depth}\t$variant[$i]->{genotype}\t$variant[$i]->{freq}";
					$samples_for_indel ++;
					$SB = $SB + $variant[$i]->{SB};
		
					# Read new line from file and obtain new variant
					$lines[$i] = readline($files[$i]);
					if (defined($lines[$i])){
						$variant[$i] = variants->create_from_tsv($lines[$i]);
					} else {
						$variant[$i] = undef;
					}
				} else {
					$line = $line . "\t-\t-\t-";
				}
			} else {
				$line = $line . "\t-\t-\t-";
			}
		}
		$line = $line ."\t" . ($SB/$samples_for_indel) . "\n";
		print OUT $line;
	}
} while (defined($first_variant));

# Close files
close OUT;

# Finish log
print LOG "\nProcess finished... ".localtime()."\n";
close LOG;
print  "\nProcess finished... ".localtime()."\n";
exit;