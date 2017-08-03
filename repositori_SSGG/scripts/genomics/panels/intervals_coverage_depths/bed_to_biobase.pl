#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: bed_to_biobase.pl
# @Author: gmarco
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use File::Path qw(make_path remove_tree);


my $biobase_path = $ENV{'BIOBASE_PATH'};
my $tmp_dir = "$biobase_path/tmp_bed";

my @biobase_files = ('bioBaseSortUniques_chr1.txt', 'bioBaseSortUniques_chr2.txt', 'bioBaseSortUniques_chr3.txt', 'bioBaseSortUniques_chr4.txt', 'bioBaseSortUniques_chr5.txt', 'bioBaseSortUniques_chr6.txt', 'bioBaseSortUniques_chr7.txt', 'bioBaseSortUniques_chr8.txt', 'bioBaseSortUniques_chr9.txt', 'bioBaseSortUniques_chr10.txt', 'bioBaseSortUniques_chr11.txt', 'bioBaseSortUniques_chr12.txt', 'bioBaseSortUniques_chr13.txt', 'bioBaseSortUniques_chr14.txt', 'bioBaseSortUniques_chr15.txt', 'bioBaseSortUniques_chr16.txt', 'bioBaseSortUniques_chr17.txt', 'bioBaseSortUniques_chr18.txt', 'bioBaseSortUniques_chr19.txt', 'bioBaseSortUniques_chr20.txt', 'bioBaseSortUniques_chr21.txt', 'bioBaseSortUniques_chr22.txt', 'bioBaseSortUniques_chrX.txt', 'bioBaseSortUniques_chrY.txt', 'bioBaseSortUniques_chrM.txt');

main();

sub main {

	my $chr;
	#my $start_pos, $end_pos, $acc_num, $gene, $dbsnp, $hgvs_cdna, $disease;
	 
	if (-d $tmp_dir){
		remove_tree($tmp_dir, {keep_root => 1} );
		sleep 1;
	}
	
	else {
		make_path($tmp_dir, { verbose => 1, mode => 0775});
		sleep 1;
	}
	
	foreach my $biobase_chr_file (@biobase_files){
		
		#open biobase chr file
		open BIOBASE_FILE, "$biobase_path/$biobase_chr_file" or die $!;
		
		#open tmp file to output bed format
		open TMP_FILE, ">$tmp_dir/tmp_$biobase_chr_file" or die $!;
	
		
		my $first = 1;
		while (<BIOBASE_FILE>){
			if ($first == 1){
				$first = 0;
				next;
			}
			
			my @line_items = split ("\t", $_);
			

			$chr = (split (":", $line_items[4]))[0];
			my $start_pos = (split (":", $line_items[4]))[1];
			my $end_pos = (split (":", $line_items[4]))[2];
			my $acc_num = $line_items[1];
			my $gene = $line_items[2];
			my $dbsnp = $line_items[6];
			my $hgvs_cdna = $line_items[9];
			my $disease = $line_items[7];
			
			#Skip malformed lines.
			if (not defined $start_pos){
				next;
			}
			
			#if position looks like: chr:start-end
			if ($start_pos =~ /-/ ){
				$end_pos = (split("-", $start_pos))[1];
				$start_pos = (split("-", $start_pos))[0];
			}
			
			#if position looks like: chr:start 
			else {
				$end_pos = $start_pos+1;
			}
			
			#chr1	123123	123124	CM071564|ATP13A2|rs121918227|NM_022089.2:c.1510G-C|Parkinsonism, juvenile|
			print TMP_FILE "$chr\t$start_pos\t$end_pos\t", join("|", $acc_num, $gene, $dbsnp, $hgvs_cdna, $disease), "\n";
			
		}
		#Sort bed file
		close TMP_FILE;
		`sort -k1,1 -k2,2g -o $tmp_dir/sorted_tmp_$biobase_chr_file $tmp_dir/tmp_$biobase_chr_file`;
		`rm $tmp_dir/tmp_$biobase_chr_file`;
		
		print "$biobase_chr_file processed !\n";
		
	}	
	
	open MASTERBIOBED, ">$biobase_path/masterbioBed_sorted.bed";
	
	print "Concatenating files into final bed..\n";
	#concatenate files into final bed 
	foreach my $biobase_tmp_chr_file (@biobase_files){
		open TMP_SORTED, "$tmp_dir/sorted_tmp_$biobase_tmp_chr_file\n" or die $!;
		while (<TMP_SORTED>){
			print MASTERBIOBED;
		}
	}
	
	
	print "Removing tmp folder..\n";
	remove_tree($tmp_dir);
	
	close MASTERBIOBED;
	chmod 0775, "$biobase_path/masterbioBed_sorted.bed";
	
	print "Successfully updated: $biobase_path/masterbioBed_sorted.bed\n";
	
}

