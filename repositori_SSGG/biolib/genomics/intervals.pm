#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Library containing intervals utilities
# @Author: Arbol
# @Contributors: biomart & ensembl scripts
########################################################

package intervals;

use strict;
use warnings;
use Data::Dumper;
use Exporter;
use Storable qw(dclone);

our @ISA = qw(Exporter);
our @EXPORT = qw(subtract_intervals_arrays add_intervals_arrays add_intervals_to_hash count_intervals_hash add_and_short_intervals subtract_intervals print_and_count_intervals add_intervals_hashes subtract_intervals_hashes normalize_intervals_hash obtain_splicing_from_exons count_intervals);

# -----------------------
# -- Include libraries --
# -----------------------


# ---------------
# -- functions --
# ---------------

sub add_intervals_arrays {
	# This function merges the intervals from 2 different intervals arrays.
	
	# Arguments:
	my ( $intervals1, $intervals2 ) = @_;
	my @output;
	
	for (my $i=0; $i <= $#{$intervals1}; $i++) {
		my @interval_to_add = ( ${$intervals1}[$i][0], ${$intervals1}[$i][1] );
		add_and_short_intervals ( \@output, \@interval_to_add );
	}
	for (my $i=0; $i <= $#{$intervals2}; $i++) {
		my @interval_to_add = ( ${$intervals2}[$i][0], ${$intervals2}[$i][1] );
		add_and_short_intervals ( \@output, \@interval_to_add );
	}
	return \@output;
}

sub subtract_intervals_arrays {
	# This function subtracts intervals2 array from intervals1 array
	
	# Arguments:
	my ( $intervals1, $intervals2 ) = @_;
	
	for (my $i=0;$i<=$#{$intervals2};$i++) {
		my @interval_to_subtract = (${$intervals2}[$i][0], ${$intervals2}[$i][1]);
		subtract_intervals ( $intervals1, \@interval_to_subtract );
	}
	return $intervals1;
}

sub show_intervals_array {
	# This is an auxiliar subroutine for showing the intervals of an array
	
	# Arguments:
	my ( $intervals_array ) = @_;
	
	#my $now = localtime;
	#print "Status at $now\n";
	
	foreach ( @{$intervals_array} ) {
		print "( " . ${$_}[0] . ", " . ${$_}[1] . " )\n";
	}
	print "\n";
}

sub obtain_splicing_from_exons {
	# This is an auxiliar subroutine for obtaining the splicing intervals array from exons (CDS) intervals arrays
	
	# Arguments:
	my ( $intervals_array, $splicing_bases ) = @_;
	my @splicing_intervals;
	
	foreach ( @{$intervals_array} ) {
		my @interval1;
		my @interval2;
		my @aux;
		if (${$_}[0] > ${$_}[1]) {
			$aux[0] = ${$_}[1];
			$aux[1] = ${$_}[0];
		} else {
			$aux[0] = ${$_}[0];
			$aux[1] = ${$_}[1];
		}
		if (${$_}[0] < $splicing_bases) {$interval1[0] = 0;}
		else {$interval1[0] = $aux[0] - $splicing_bases;}
		$interval1[1] = $aux[0] - 1;
		$interval2[0] = $aux[1] + 1;
		$interval2[1] = $aux[1] + $splicing_bases;
		add_and_short_intervals(\@splicing_intervals,\@interval1);
		add_and_short_intervals(\@splicing_intervals,\@interval2);
	}
	return \@splicing_intervals;
}

sub add_and_short_intervals {
	# This subroutine shorts and merges new intervals into an existent intervals array.
	# Note:
	#          low2  high2
	# new  :     [     ]
	#          low1  high1
	# array:           0           1                               2                               3                
	#             [     ]  [               ]            [                   ]               [              ]        
	# Cases of new interval:
	# A                           0   0
	#                             [   ]
	#                             1   1
	# B                           0             1
	#                             [             ]
	#                             1             1
	# C                           0                           0
	#                             [                           ]
	#                             1                           2
	# D                           0                                                   1
	#                             [                                                   ]
	#                             1                                                   2
	# E                                         1    1
	#                                           [    ]
	#                                           2    1
	# F                                         1                0
	#                                           [                ]
	#                                           2                2
	# G                                         1                                       1
	#                                           [                                       ]
	#                                           2                                       2
	# H                                         1                                                   0
	#                                           [                                                   ]
	#                                           2                                                   3
	# I 1   
	#   [    ]
	#   0   -1
	# J                                                                                                                     
	#                                                                                                               [       ]
	#                                                                                                               4       
	
	# Arguments:
	# Please note that $new_interval, $intervals_array and $added_pos are all references!!
	my ( $intervals_array, $new_interval ) = @_;
	
	# Normalize current interval (low,high):
	my @current_interval;
	#unless ( looks_like_number( ${$new_interval}[0]) ) { print "${$new_interval}[0] \n";}
	if ( ${$new_interval}[1] > ${$new_interval}[0] ){ @current_interval = (${$new_interval}[0], ${$new_interval}[1]); }
	else { @current_interval = (${$new_interval}[1], ${$new_interval}[0]); }
							
	# The exon coordinates intervals get ordered in the @exon_intervals array.
	# These intervals get also merged with adjacent intervals
						
	# The first interval is directly added to the array
	if ( $#{$intervals_array} == -1 ) {
		push (@{$intervals_array}, \@current_interval);
	} else {
		# Check if the exon interval intersects with existent intervals,
		# and in which position it should be added

		# The needed variables are set (see schema below): 
		my $low_pos = 0;
		my $low_out = 0;
		while ( $low_pos <= $#{$intervals_array} && $current_interval[0] > ${$intervals_array}[$low_pos][1]+1 ) { $low_pos ++;	}
		# Case J:
		if ( $low_pos > $#{$intervals_array} ) { push(@{$intervals_array}, \@current_interval); }
		else {
			if ( $current_interval[0] < ${$intervals_array}[$low_pos][0] ) {
				$low_out = 1;
			}
			
			my $high_pos = $#{$intervals_array};
			my $high_out = 0;
			while ( $high_pos >= 0 && $current_interval[1]+1 < ${$intervals_array}[$high_pos][0] ) { $high_pos --; }
			# Case I:
			if ( $high_pos < 0 ) { unshift(@{$intervals_array}, \@current_interval); }
			else {
				if ( $current_interval[1] > ${$intervals_array}[$high_pos][1] ) {
					$high_out = 1;
				}
				# New interval doesn't intersects existent intervals, or it engulfs one or more intervals completely:
				if ($low_out and $high_out){
					splice(@{$intervals_array}, $low_pos, ($high_pos-$low_pos+1), \@current_interval);
				}
				# Both new interval boundaries are within existent intervals 
				elsif ($low_out == 0 and $high_out == 0){
					${$intervals_array}[$low_pos][1] = ${$intervals_array}[$high_pos][1];
					splice(@{$intervals_array}, $low_pos+1, ($high_pos-$low_pos));
				} else {
					# Low position of new interval is out of existent intervals
					if ($low_out){
						${$intervals_array}[$high_pos][0] = $current_interval[0];
						# Delete intervals that have been engulfed by new interval:
						splice(@{$intervals_array}, $low_pos, ($high_pos-$low_pos));
					}
					# High position of new interval is out of existent intervals
					if ($high_out){
						${$intervals_array}[$low_pos][1] = $current_interval[1];
						# Delete intervals that have been engulfed by new interval:
						splice(@{$intervals_array}, $low_pos+1, ($high_pos-$low_pos));
					}	
				}
			}
		}
	}
}

sub subtract_intervals {
	# This subroutine subtracts an interval from an existent intervals array.
	# Note: see the explanation at the beginning of add_and_short_intervals
	
	# Arguments:
	my ( $intervals_array, $subtract_interval ) = @_;
	
	# Normalize current interval (low,high):
	my @current_interval;
	if ( ${$subtract_interval}[1] > ${$subtract_interval}[0] ){ @current_interval = (${$subtract_interval}[0], ${$subtract_interval}[1]); }
	else { @current_interval = (${$subtract_interval}[1], ${$subtract_interval}[0]); }
	
	#print "current: $current_interval[0] $current_interval[1] \n";
	
	# If no interval is present in the array, nothing will be done
	if ( $#{$intervals_array} != -1 ) {
		# Check if the exon interval intersects with existent intervals,
		# and in which position it should be subtracted

		# The needed variables are set (see schema below, in add_and_short_intervals): 
		my $low1 = 0;
		my $low2;
		while ( $low1 <= $#{$intervals_array} && $current_interval[0] > ${$intervals_array}[$low1][1] ) { $low1 ++;	}
		# All cases but J:
		if ( $low1 <= $#{$intervals_array} ) {
			if ( $current_interval[0] > ${$intervals_array}[$low1][0] ) {
				$low2 = $low1 + 1;
			} else { $low2 = $low1; }
			
			my $high1 = $#{$intervals_array};
			my $high2;
			while ( $high1 >= 0 && $current_interval[1] < ${$intervals_array}[$high1][0] ) {$high1 --;}
			#print "$high1 >= 0 && $current_interval[1] < ${$intervals_array}[$high1][0] \n";
			# All cases but I:
			if ( $high1 >= 0 ) {
				if ( $current_interval[1] < ${$intervals_array}[$high1][1] ) {
					$high2 = $high1 - 1;
				} else { $high2= $high1; }
				# Cases E and G:
				if ( $low1 == $low2 && $high1 == $high2 ) {
					splice(@{$intervals_array}, $low1, ($high1-$low1+1));
				}
				else {
					# Cases F and H:
					if ( $low1 == $low2 ) {
						${$intervals_array}[$high1][0] = $current_interval[1]+1;
					# Cases B and D:
					} elsif ( $high1 == $high2 ) {
						${$intervals_array}[$low1][1] = $current_interval[0]-1;
					} 
					# Cases C, D and H:
					if ( $high1 > $low1 ) {
						# Cases C and D:
						if ( $low1 != $low2 ) { 
							${$intervals_array}[$low1][1] = $current_interval[0]-1;
						}
						# Cases C and H:
						if ($high1 != $high2) {
							#print "C y H\n";
							${$intervals_array}[$high1][0] = $current_interval[1]+1;
						}
						# Cases C,D and H:
						splice(@{$intervals_array}, $low2, ($high2-$low2+1));
					}
					# Case A:
					if ($high2<$low1 && $low1!=$low2) {
						my @new_interval = ($current_interval[1]+1, ${$intervals_array}[$low1][1]);
						${$intervals_array}[$low1][1] = $current_interval[0]-1;
						splice(@{$intervals_array}, $low2, 0, \@new_interval);
					}
				}
			}
		}
	}
}

sub count_intervals {
	# This function counts the positions belonging to a list of intervals
	
	# Arguments:
	my ( $intervals ) = @_;
	
	# Vars definitions:
	my $sum_total = 0;
	
	# Count intervals from array:
	for (my $i=0; $i < $#{$intervals} + 1; $i++) {
		if (length(${$intervals}[$i]) > 0){ 
			$sum_total = $sum_total + ${$intervals}[$i][1] - ${$intervals}[$i][0] + 1;
		}
	}
	return $sum_total;
}

sub print_and_count_intervals {
	# This function counts the positions belonging to a list of intervals and outputs its information to a file
	
	# Arguments:
	my ( $intervals, $chromosome, $file, $append_mode ) = @_;
	
	# Vars definitions:
	my $sum_total = 0;
	
	# File output
	if ( $append_mode ){
		open (FILE, ">>", $file) or die("ERROR: Couldn't open file $file in append mode!! $!\n");
	}
	else{
		open (FILE, ">", $file) or die("ERROR: Couldn't open output file $file!! $!\n");
	}
	
	# Count and print intervals from array:
	if (ref($intervals) eq 'ARRAY'){
		for (my $i=0; $i <= $#{$intervals}; $i++) {
			my $out_row = "chr" . $chromosome . "\t" . ${$intervals}[$i][0] . "\t" . ${$intervals}[$i][1] . "\n";
			print FILE $out_row;
			$sum_total = $sum_total + ${$intervals}[$i][1] - ${$intervals}[$i][0] + 1;
		}
		close FILE;
		return $sum_total;
	} else {
		return 0;
	}
}

sub show_intervals_hash {
	# This is an auxiliar subroutine for showing any kind of hash
	
	# Arguments:
	my ( $hash ) = @_;
	
	foreach (keys %$hash){
		my $one_hash_key = $_;
		print "key: $one_hash_key\n"; 
		
		if (ref($$hash{$one_hash_key}) eq "HASH"){
			show_intervals_hash ($hash->{$one_hash_key});
		}
		else {
			show_intervals_array ( \@{$hash->{$one_hash_key}} );
		}
	}
}

sub add_intervals_to_hash{
	# This function adds elements to hashes, with as much key-levels as desired and checking if the keys already exist
	
	# Arguments:
	my ( $hash, $intervals, @hash_keys ) = @_;
	
	if ( $#hash_keys < 0 ){
		die ("ERROR: No hash keys passed as argument. Impossible to add to hash without a given key!!\n");
	}
	
	# Search the key in the hash:
	my $found = 0;
	foreach ( keys %$hash ) {
		my $one_hash_key = $_;
		# Add exon intervals:
		if ( $hash_keys[0] eq $one_hash_key ){
			if ($#hash_keys > 0){
				$hash->{$one_hash_key} = add_intervals_to_hash (\%{$$hash{$one_hash_key}}, \@{$intervals}, @hash_keys[1..$#hash_keys]);
			} else{
				$hash->{$one_hash_key} = add_intervals_arrays (\@{$$hash{$one_hash_key}}, \@{$intervals});
			}
			$found = 1;
		}
	}
	# If the hash key wasn't found in the hash:
	if ( $found == 0 ) {
		if ($#hash_keys > 0){
			$hash->{$hash_keys[0]} = add_intervals_to_hash (\%{$$hash{$hash_keys[0]}}, \@{$intervals}, @hash_keys[1..$#hash_keys]);
		} else{
			$hash->{$hash_keys[0]} = \@$intervals;
		}
	}
	return $hash;
}

sub add_intervals_hashes {
	# This function merges the intervals from 2 different intervals hashes.
	
	# Arguments:
	my ( $intervals1, $intervals2 ) = @_;
	
	# Vars definitions:
	my %output = ();
	
	if ( keys(%$intervals1) > 0 ) {
		foreach ( keys %$intervals1 ) {
			my $stored_chrom = $_;
			if (ref($intervals1->{$stored_chrom}) eq "HASH"){
				foreach ( keys %{$intervals1->{$stored_chrom}}){
					add_intervals_to_hash (\%output, \@{$intervals1->{$stored_chrom}->{$_}}, $stored_chrom, $_);
				}
			} else {
				add_intervals_to_hash (\%output, \@{$intervals1->{$stored_chrom}}, $stored_chrom);
			}
		}
	}
	if ( keys(%$intervals2) > 0 ) {
		foreach ( keys %$intervals2 ) {
			my $stored_chrom = $_;
			if (ref($intervals2->{$stored_chrom}) eq "HASH"){
				foreach ( keys %{$intervals2->{$stored_chrom}}){
					add_intervals_to_hash (\%output, \@{$intervals2->{$stored_chrom}->{$_}}, $stored_chrom, $_);
				}
			} else {
				add_intervals_to_hash (\%output, \@{$$intervals2{$stored_chrom}}, $stored_chrom);
			}
		}
	}
	return \%output;
}

sub subtract_intervals_from_hash
# This function subtracts elements from hashes, with as much key-levels as desired and checking if the keys already exist
{
	# Arguments:
	my ($hash, $intervals, @hash_keys) = @_;
	
	if ( $#hash_keys < 0 ){
		die ("ERROR: No hash keys passed as argument. Impossible to add to hash without a given key!!\n");
	}
	
	# Search the key in the hash:
	foreach my $one_hash_key (keys %$hash){
		# Add exon intervals:
		if ( $hash_keys[0] eq $one_hash_key ){
			if ($#hash_keys > 0){
				$hash->{$one_hash_key} = subtract_intervals_from_hash (\%{$$hash{$one_hash_key}}, \@{$intervals}, @hash_keys[1..$#hash_keys]);
			} else{
				if (ref($hash->{$one_hash_key}) eq "HASH"){
					foreach (keys %{$hash->{$one_hash_key}}){
						$hash->{$one_hash_key} = subtract_intervals_from_hash (\%{$hash->{$one_hash_key}}, \@{$intervals}, $_);
					}
				}
				else {
					subtract_intervals_arrays (\@{$hash->{$one_hash_key}}, \@{$intervals});
				}
			}
		}
	}
	# If the hash key wasn't found in the hash then nothing will be subtracted
	return $hash;
}

sub subtract_intervals_hashes {
	# This function subtracts the intervals from a hash (intervals2) to another (intervals1) and returns the non-overlapping intervals
	
	# Arguments:
	my ( $intervals1, $intervals2 ) = @_;
	
	# Vars definitions:
	# $intervals1 must be copied before starting
	my $output = dclone($intervals1);
	
	if ( keys(%$intervals2) > 0 ) {
		foreach my $stored_chrom (keys %$intervals2) {
			subtract_intervals_from_hash ($output,\@{$$intervals2{$stored_chrom}},$stored_chrom);
			if (defined($$output{$stored_chrom})){
				if (ref($$output{$stored_chrom}) eq "HASH"){
					foreach ( keys %{$$output{$stored_chrom}}){
						if (not defined($$output{$stored_chrom}{$_}[0])){
							delete($$output{$stored_chrom}{$_});
						}
					}
				}
				if ($$output{$stored_chrom} eq ""){
					delete($$output{$stored_chrom});
				}
			}
		}
	}
	return $output;
}

sub normalize_interval {
# This function returns a normalized interval (begin coordinate smaller than final coordinate)
	# Arguments:
	my ($interval) = @_;
	
	# Normalize interval:
	if (${$interval}[0] > ${$interval}[1]){
		my $aux = ${$interval}[0];
		${$interval}[0] = ${$interval}[1];
		${$interval}[1] = $aux;
	}
	return $interval;
}

sub intersect_intervals {
# This function returns the percentage and amount of bases of the intersection between intervals
	
	# Arguments:
	my ($interval1, $interval2, $intersection) = @_;
	
	# Vars definition:
	my @result = (0,0); # Results vector, containing percentage and amount of bases of the intersection
	
	# Normalize intervals:
	$interval1 = normalize_interval($interval1);
	$interval2 = normalize_interval($interval2);
	
	# Check intersection:
	if (${$interval1}[0] <= ${$interval2}[0] and ${$interval1}[1] >= ${$interval2}[0]){
		if (${$interval1}[1] < ${$interval2}[1]){
			$result[1] = ${$interval1}[1] - ${$interval2}[0] + 1;
		} else {
			$result[1] = ${$interval2}[1] - ${$interval2}[0] + 1;
		}
		$result[0] = $result[1]/(${$interval1}[1] - ${$interval1}[0] + 1);
	}
	if (${$interval1}[0] >= ${$interval2}[0] and ${$interval1}[0] <= ${$interval2}[1]){
		if (${$interval2}[1] < ${$interval1}[1]){
			$result[1] = ${$interval2}[1] - ${$interval1}[0] + 1;
		} else {
			$result[1] = ${$interval1}[1] - ${$interval1}[0] + 1;
		}
		$result[0] = $result[1]/(${$interval1}[1] - ${$interval1}[0] + 1);
	}
	return @result;
}

sub filtering_repeats {
# This function filters out the intervals from a hash that overlap with other intervals in the same hash
	
	# Arguments:
	my ( $hash,$bases_threshold,$LOG ) = @_;
	
	# Vars definition:
	my @result;
	
	# Take each interval from hash1:
	my @keys_hash = keys %$hash;

	# get all the values, in order
	my $index1 = 0;
	while($index1 <= $#keys_hash){
		my $one_hash_key = $keys_hash[$index1];
		my @keys_miRNA = keys %{$hash->{$one_hash_key}};
		my $chr1 = $keys_miRNA[0];
		my $interval1 = ${$hash->{$one_hash_key}->{$chr1}}[0];
		my $index2 = $index1 + 1;
		my $intersection = 0;
		while($index2 <= $#keys_hash){
			my $other_hash_key = $keys_hash[$index2];
			if ($one_hash_key ne $other_hash_key){
				my @keys_miRNA2 = keys %{$hash->{$other_hash_key}};
				my $chr2 = $keys_miRNA2[0];
				if ($#{$hash->{$other_hash_key}{$chr2}} >= 0){
					my $interval2 = ${$hash->{$other_hash_key}->{$chr2}}[0];
					#print "chr1 $chr1 int1 $interval1 : (${$interval1}[0] , ${$interval1}[1]) chr2 $chr2 int2 $interval2 : (${$interval2}[0] , ${$interval2}[1])\n";
					if ($chr1 eq $chr2){
						@result = intersect_intervals ($interval1, $interval2);	
					}
					else {
						@result = (0,0)
					}
					if ($result[1] >= $bases_threshold){
						my $percentage = $result[0] * 100; 
						my $line = "INFO: Interval ( ${$interval1}[0], ${$interval1}[1] ) intersects interval ( ${$interval2}[0], ${$interval2}[1] ) with $percentage % and $result[1] bp. Intervals get merged and counts summed up.\n";
						print $line;
						print $LOG $line;
						add_intervals_to_hash($hash->{$other_hash_key},$hash->{$one_hash_key}->{$chr1},$chr1);
						my @values1 = split (';',$one_hash_key);
						my @values2 = split (';',$other_hash_key);
						my @values;
						if ($values1[1] > $values2[1]){
							$values[0] = $values1[0];
							$values[1] = $values1[1];
							$values[4] = $values1[4];
							$values[5] = $values1[5];
						} else {
							$values[0] = $values2[0];
							$values[1] = $values2[1];
							$values[4] = $values2[4];
							$values[5] = $values2[5];
						}
						$values[2] = $values1[2] + $values2[2];
						$values[3] = abs($hash->{$other_hash_key}{$chr1}[0][1]-$hash->{$other_hash_key}{$chr1}[0][0]);
						
						my $new_key = $values[0] . ";". $values[1] . ";" . $values[2] . ";" . $values[3] . ";" . $values[4] . ";" . $values[5];
						
						$hash->{$new_key} = delete $hash->{$other_hash_key};
						splice(@keys_hash,$index2,1,$new_key);
						$intersection = 1;
					}
				}
			}
			$index2++;
		}
		if ($intersection){
			delete $hash->{$one_hash_key};
		}
		else {
			$index1++;	
		}
	}
}

sub overlapping_intervals_hashes {
# This function detects which intervals in hash1 overlap intervals in hash2, and dismish them from hash1 according to the defined criteria
# $comparison is the kind of comparison between the intervals intersection and the thresholds. Its values are: le (lesser or equal), ge (greater or equal)
	
	# Arguments:
	my ( $hash1,$hash2,$percentage_threshold,$bases_threshold,$comparison,$LOG ) = @_;
	
	# Vars definition:
	my @result;
	
	# Take each interval from hash1:
	foreach (keys %$hash1){
		my $one_hash_key = $_;
		if (ref($$hash1{$one_hash_key}) eq "HASH"){
			overlapping_intervals_hashes ($hash1->{$one_hash_key},$hash2,$percentage_threshold,$bases_threshold,$comparison,$LOG);
		}
		else {
			my $i1 = 0;
			while ($i1<=$#{$hash1->{$one_hash_key}}){
				my $interval1 = ${$hash1->{$one_hash_key}}[$i1];
				my $increment = 1;
				#Take each interval from hash2:
				foreach (keys %$hash2){
					my $another_hash_key = $_;
					if (ref($$hash2{$another_hash_key}) eq "HASH"){
						overlapping_intervals_hashes ($hash1,$hash2->{$another_hash_key},$percentage_threshold,$bases_threshold,$comparison,$LOG);
					}
					# Intervals from same key (usually chromosome) will be checked
					elsif ($one_hash_key eq $another_hash_key){
						foreach (@{$hash2->{$another_hash_key}}){
							my $interval2 = $_;
							$increment = 0;
							@result = intersect_intervals ($interval1, $interval2);
							my $dismish = 0;
							if ($comparison eq "le"){
								if ($result[0] <= $percentage_threshold or $result[1] <= $bases_threshold){
									$dismish = 1;
								}
								else {
									$increment = 1;
								}
							} elsif ($comparison eq "ge"){
								if ($result[0] >= $percentage_threshold or $result[1] >= $bases_threshold){
									$dismish = 1;
								}
								else {
									$increment = 1;
								}
							} else {$increment = 1;}
							if ($dismish){
								my $percentage = $result[0] * 100; 
								my $line = "INFO: Interval ( ${$interval1}[0], ${$interval1}[1] ) intersects interval ( ${$interval2}[0], ${$interval2}[1] ) with $percentage % and $result[1] bp\n";
								print $line;
								print $LOG $line;
								splice(@{$hash1->{$one_hash_key}}, $i1, 1);
							}
						}
					}
				}
				if ($increment){
					$i1++;
					$increment = 0;
				}
			}
		}
	}
}

sub count_intervals_hash {
	# This function counts the positions belonging to an intervals hash
	
	# Arguments:
	my ( $intervals ) = @_;
	
	# Vars definitions:
	my @chroms = (1..22, "X", "Y", "M");
	my $sum_total = 0;
	
	foreach ( @chroms ) {
		my $chrom = $_;
		if ( exists($$intervals{$chrom}) ) {
			if (ref($intervals->{$chrom}) eq "HASH"){
				foreach ( keys %{$intervals->{$chrom}}){
					$sum_total = $sum_total + count_intervals(\@{$intervals->{$chrom}->{$_}});
				}
			} else {
				$sum_total = $sum_total + count_intervals($$intervals{$chrom});
			}
		}
	}
	return $sum_total;
}

sub print_and_count_intervals_hash {
	# This function counts the positions belonging to an intervals hash and outputs its information to a file
	
	# Arguments:
	my ( $intervals, $file ) = @_;
	
	# Vars definitions:
	my @chroms = (1..22, "X", "Y", "M");
	my $first = 0;
	my $sum_total = 0;
	
	foreach ( @chroms ) {
		my $chrom = $_;
		if ( exists($$intervals{$chrom}) ) {
			if ($first==0){$sum_total = print_and_count_intervals(@{$$intervals{$_}},$chrom,$file,0);}
			else {$sum_total = print_and_count_intervals(@{$$intervals{$_}},$chrom,$file,1) + $sum_total;}
		}
	}
	return $sum_total;
}

sub normalize_intervals_hash {
	# This subroutine normalizes an intervals hash, so that its intervals are at least equal or bigger than $minsize
	
	# Arguments:
	my ( $intervals, $LOG, $minsize ) = @_;
	
	# Parse $intervals
	my @chroms = (1..22, "X", "Y", "M");
	foreach ( @chroms ) {
		my $current_chrom = $_;
		if ( exists($$intervals{$_}) ) {
			for (my $i=0;$i<=$#{$$intervals{$current_chrom}};$i++){
				# Normalization of all intervals for eArray (Agilent's utility for designing panels)
				# Exon intervals under $minsize (default 120nt) will be "enlarged" until reaching this size
				my $intervals_difference = $minsize - ${$$intervals{$current_chrom}}[$i][1] + ${$$intervals{$current_chrom}}[$i][0];
				if ( $intervals_difference > 0 ){
					my $suplement = int ( $intervals_difference / 2); # + $intervals_difference % 2;
					my @new_interval;
					$new_interval[0] = ${$$intervals{$current_chrom}}[$i][0] - $suplement;
					$new_interval[1] = ${$$intervals{$current_chrom}}[$i][1] + $suplement;
					my $line = "INFO: Interval (chr" . $current_chrom . ":" . ${$$intervals{$current_chrom}}[$i][0] . "," . ${$$intervals{$current_chrom}}[$i][1] . ") smaller than $minsize!! It will be normalized to (" . $new_interval[0] . "," . $new_interval[1] . ")!!\n";
					print $line;
					print $LOG $line;
					add_and_short_intervals ( \@{$$intervals{$current_chrom}}, \@new_interval );
				}
			}
		}
	}
}

1;