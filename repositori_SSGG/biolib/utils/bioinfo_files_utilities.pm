#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Library containing file import & export utilities
# @Author: Arbol
# @Contributors: biomart & ensembl scripts
########################################################

package bioinfo_files_utilities;

use strict;
use warnings;
use Data::Dumper;
use Exporter;
use intervals qw(add_intervals_to_hash count_intervals_hash);

our @ISA = qw(Exporter);
our @EXPORT = qw(import_intervals_file create_intervals_file);

# ---------------
# -- functions --
# ---------------

sub import_intervals_file {
# This is an auxiliar subroutine for importing bed, gff, gff3, tsv and chr-coord files into an intervals hash
	
	# Arguments:
	my ( $file, $file_type, $intervals, $LOG, $get_name) = @_;
	
	# File open:
	open (FILE, "<", $file) or die("ERROR: Couldn't open input file $file!! $!\n");
	
	# Variables:
	my $gff3_chr = "";
	
	# Count and print intervals from array:
	while ( <FILE> ){
		chomp($_);
		if ($file_type eq "bed"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[1]>$values[2]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				
				# Get also interval name (in case that variant name is included)
				my @interval = [($values[1],$values[2])];
				if (defined($get_name) and defined($values[3])){
					add_intervals_to_hash($intervals, \@interval, $tmp[1], $values[3]);
				} else {
					add_intervals_to_hash($intervals, \@interval, $tmp[1]);
				}
			}
		} elsif ($file_type eq "gff"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[3]>$values[4]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				my @interval = [($values[3],$values[4])];
				add_intervals_to_hash($intervals, \@interval, $tmp[1]);
			}
		} elsif ($file_type eq "gff3"){
			my $file_line = $_;
			if (substr($file_line,0,2) ne "##"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[3]>$values[4]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				my @interval = [($values[3],$values[4])];
				if ($gff3_chr eq ""){
					my $line = "ERROR: Chromosome not defined. Interval cannot be stored into hash!!\n";
					print $line;
					print $LOG $line;
				} else {
					add_intervals_to_hash($intervals, \@interval, $gff3_chr);
				}
			}
			else {
				my $chrom_string = "chromosome";
				if (index($file_line, $chrom_string) != -1){
					my @tmp2 = split (' ', $file_line);
					my @tmp3 = split (',', $tmp2[4]);
					$gff3_chr = $tmp3[0];
				}
			}
		} elsif ($file_type eq "tsv"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[1]>$values[2]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				my @interval = [($values[1],$values[2])];
				add_intervals_to_hash($intervals, \@interval, $tmp[1]);
			}
		} elsif ($file_type eq "chr-coord"){
			my @values = split (':', $_);
			my @tmp = split('chr',$values[0]);
			my @interval = [split('-',$values[1])];
			if (defined($values[2])){
				add_intervals_to_hash($intervals,\@interval,$tmp[1],$values[2]);
			} else {
				add_intervals_to_hash($intervals, \@interval, $tmp[1]);
			}
		}
	}
	close FILE;

	# Count the bases in the created intervals hash:
	return count_intervals_hash($intervals);
}

sub output_to_file {
# This is an auxiliar subroutine for output data to files with formats bed, gff, gff3, tsv and chr-coord
	
	# Arguments:
	my ( $FILE, $file_type, $intervals, $LOG) = @_;
	
	# File open:
	open (FILE, "<", $FILE) or die("ERROR: Couldn't open input file $FILE!! $!\n");
	
	# Variables:
	my $gff3_chr = "";	
}

sub file_conversion {
# This is an auxiliar subroutine for converting files from different formats (bed, gff, gff3, tsv and chr-coord)
	
	# Arguments:
	my ( $file1, $file1_type, $file2, $file2_type, $LOG) = @_;
	
	# File open:
	open (FILE1, "<", $file1) or die("ERROR: Couldn't open input file $file1!! $!\n");
	open (my $FILE2, ">", $file2) or die("ERROR: Couldn't open output file $file2!! $!\n");
	
	# Variables:
	my $gff3_chr = "";
	
	# Count and print intervals from array:
	while ( <FILE1> ){
		chomp($_);
		if ($file1_type eq "bed"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[1]>$values[2]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				output_to_file($FILE2,$file2_type,$tmp[1],$values[1],$values[2],$LOG);
			}
		} elsif ($file1_type eq "gff"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[3]>$values[4]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				output_to_file($FILE2,$file2_type,$tmp[1],$values[3],$values[4],$LOG);
			}
		} elsif ($file1_type eq "gff3"){
			my $file_line = $_;
			if (substr($file_line,0,2) ne "##"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[3]>$values[4]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				my @interval = [($values[3],$values[4])];
				if ($gff3_chr eq ""){
					my $line = "ERROR: Chromosome not defined. Interval cannot be stored into hash!!\n";
					print $line;
					print $LOG $line;
				} else {
					output_to_file($FILE2,$file2_type,$gff3_chr,$values[3],$values[4],$LOG);
				}
			}
			else {
				my $chrom_string = "chromosome";
				if (index($file_line, $chrom_string) != -1){
					my @tmp2 = split (' ', $file_line);
					my @tmp3 = split (',', $tmp2[4]);
					$gff3_chr = $tmp3[0];
				}
			}
		} elsif ($file1_type eq "tsv"){
			if (substr($_,0,3) eq "chr"){
				my @values = split ('\t', $_);
				my @tmp = split('chr',$values[0]);
				if ($values[1]>$values[2]){
					my $line = "WARNING: Abnormal interval in eArray desing!! Initial coordinate > End coordinate!!\n";
					print $line;
					print $LOG $line;
				}
				output_to_file($FILE2,$file2_type,$tmp[1],$values[1],$values[2],$LOG);
			}
		} elsif ($file1_type eq "chr-coord"){
			my @values = split (':', $_);
			my @tmp = split('chr',$values[0]);
			my @interval = [split('-',$values[1])];
			if (defined($values[2])){
				output_to_file($FILE2,$file2_type,$tmp[1],$interval[0],$interval[1],$LOG);
			} else {
				output_to_file($FILE2,$file2_type,$tmp[1],$interval[0],$interval[1],$LOG);
			}
		}
	}
	close FILE1;
	close $FILE2;
}


sub import_miRDeep2_file {
	# This is an auxiliar subroutine for importing miRDeep miRNAs files into an intervals hash
	
	# Arguments:
	my ( $file, $intervals, $LOG) = @_;
	
	# File open:
	open (FILE, "<", $file) or die("ERROR: Couldn't open input file $file!! $!\n");
	
	# Count and print intervals from array:
	while ( <FILE> ){
		chomp($_);
		my $line = $_;
		my $tmp1 = substr($line,0,3);
		if ($tmp1 eq "chr"){
			my @values = split ('\t', $line);
			my @tmp2 = split('chr',$values[16]);
			my @tmp3 = split(':',$tmp2[1]);
			my @interval = [split('\.\.',$tmp3[1])];
			my $interval_size = abs($interval[0][1] - $interval[0][0]);
			add_intervals_to_hash($intervals,\@interval,"$values[0];$values[1];$values[4];$interval_size;$tmp3[2];$values[15]",$tmp3[0]);
		}
	}
	close FILE;
}

sub create_intervals_file {
# This subroutine parses an intervals hash and prints it to a file with following formats: chr-coord, gff, bed
	
	# Arguments:
	my ( $output_file, $format, $LOG, $intervals, $args, $append) = @_;
	
	# Vars definition:
	my $count = 0;
	
	# Open file for writing or appending
	if ($append){
		open SEL, '>>', $output_file or die ("ERROR: Impossible to append to $output_file!! $!");
	} else {
		open SEL, '>', $output_file or die ("ERROR: Impossible to create $output_file!! $!");
	}
	
	# Print Header for gff format
	if ($format eq "gff"){
		print SEL "Chromosome,StartCoordinate,StopCoordinate,TargetType,Density,Labels\n";
	}
	
	# Check if the hash is in the chromosomes level:
	my @chroms = (1..22, "X", "Y", "M");
	my $found = 0;
	foreach ( @chroms ) {
		if (exists($intervals->{$_})){
			$found = 1;
		}
	}
	my $new_args;
	if ($found == 0){
		# If the hash is a hash of hashes with the form: hash{key1}{key2}...{keyN}{chr}
		foreach (keys %$intervals){
			my $one_hash_key = $_;
			if (ref($intervals->{$one_hash_key}) eq "HASH"){
				if ($args ne ""){
					$new_args = $args . ";" . $one_hash_key;
				} else {
					$new_args = $one_hash_key;
				}
				create_intervals_file($output_file,$format,$LOG,$intervals->{$one_hash_key},$new_args,1);
			} else {
				if ($one_hash_key ne ""){
					my $line = "WARNING: Apparently not an intervals hash ($one_hash_key)!!\n";
					print $line;
					print $LOG $line;
				}
			}
		}
	} else {
		# Now it's an intervals hash, and it can be parsed:
		foreach ( @chroms ) {
			my $current_chrom = $_;
			if ( exists($$intervals{$current_chrom})) {
				if (ref($intervals->{$current_chrom}) eq "HASH"){
					foreach (keys %{$intervals->{$current_chrom}}){
						my $one_hash_key=$_;
						foreach ( @{$intervals->{$current_chrom}->{$one_hash_key}} ){
							if ($format eq "chr-coord"){
								# output to file
								my $line = "chr" . $current_chrom . ":" . ${$_}[0] . "-" . ${$_}[1] . ":" . $one_hash_key . "\n";
								print SEL $line;
							} elsif ($format eq "gff"){
								# output to file
								my $line = $current_chrom . "," . ${$_}[0] . "," . ${$_}[1] . ",FullRegion,Standard,SSGG_Panel," . $one_hash_key . "\n";
								print SEL $line;					
							} elsif ($format eq "bed"){
								# output to file
								my $line = "chr" . $current_chrom . "\t" . ${$_}[0] . "\t" . ${$_}[1] . "\t" . $one_hash_key . "\n";
								print SEL $line;
							}
							# Count bases in intervals hash:
							$count = $count + ${$_}[1] - ${$_}[0] + 1;
						}
					}
				}
				else{
					foreach my $interval ( @{$$intervals{$current_chrom}} ){
						if ($format eq "chr-coord"){
							# output to file
							my $line = "chr" . $current_chrom . ":" . $$interval[0] . "-" . $$interval[1];
							if (not defined($args)){
								$line = $line . "\n";
							} else {
								$line = $line . ":" . $args . "\n";
							}
							print SEL $line;
						} elsif ($format eq "gff"){
							# output to file
							my $line = $current_chrom . "," . $$interval[0] . "," . $$interval[1] . ",FullRegion,Standard,SSGG_Panel";
							if (not defined($args)){
								$line = $line . "\n";
							} else {
								$line = $line . $args . "\n";
							}
							print SEL $line;
						} elsif ($format eq "bed"){
							# output to file
							my $line = "chr" . $current_chrom . "\t" . $$interval[0] . "\t" . $$interval[1];
							if (not defined($args)){
								$line = $line . "\n";
							} else {
								$line = $line . "\t" . $args . "\n";
							}
							print SEL $line;
						}
						# Count bases in intervals hash:
						$count = $count + $$interval[1] - $$interval[0] + 1;
					}
				}
			} 
		}
	}
	
	close SEL;
	return $count;
}

# **************************************************************************************************
# **************************************************************************************************
1; 