#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Copy;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

sub get_chr ($$$);
sub check_end ($$$);

my $input_id;
my $seq_id;
my $output_id;

GetOptions (
				"i=s"	=> \$input_id,	# input file
				"o=s"	=> \$output_id, # output name			
);

print "\n\n=================================================================
Parsing SAM file from BWA v 1.0\n";

if (!$input_id || !$output_id) {
	die "
		Options:
	-i BWA SAM file name
	-o Output name root 
\n".localtime()."\n=================================================================\n\n";
}

print "\nStarting process... ".localtime()."\n\n";

open (SAM, "< $input_id") or die "Can't open $input_id\n";

open (OUT, "> $output_id\_filtered.sam") or die "Can't create $output_id\_filtered.sam\n";
open (DIS, "> $output_id\_discarded.sam") or die "Can't create $output_id\_discarded.sam\n";

my %Chrm;
my %bad_id;

while (my $line = <SAM>) {
	chomp $line;
	if ($line =~ m/^\@HD/) { next;}
	elsif ($line =~ m/^\@SQ/) {
		print OUT "$line\n";
		print DIS "$line\n";
		my ($chr, $end);
		get_chr($line, \$chr, \$end);
		$Chrm{$chr}{end} = $end;
		#print "$chr\t$end\t$Chrm{$chr}{end}\n";
	}
	elsif ($line =~ m/^\@PG/ || $line =~ m/^\@RG/) { 
		print OUT "$line\n";
		print DIS "$line\n";
	}
	else {
		my @data = split (/\s+/, $line);
		my $chr = $data[2];
		my $pos = $data[3];
		if ($pos == 0 && $data[6] ne '*') {
			$data[3] = 1;
		}
		$line = join("\t", @data);
		my $length = length($data[9]);
		my $cygar = $data[5];
		my $final = &get_length($pos, $length, $cygar);
		my $pair_pos = $data[7];
		my $pair_chr;
		if ($data[6] ne '*') {
			if ($data[6] eq '=') {
				$pair_chr = $chr;
			}
			else {
				$pair_chr = $data[6];
			}
			if ($bad_id{$data[0]}) {
				print DIS "$line\n";
			}
			elsif ($chr eq '*') {
				next;
			}
			else {
				if (&check_end($chr, $final, \%Chrm) == 1 && &check_end($pair_chr, $pair_pos, \%Chrm) == 1) {
					print OUT "$line\n";
				}
				else {
					print DIS "$line\n";
					$bad_id{$data[0]} = 1;
				}
			}
		}
		else {
			if ($bad_id{$data[0]}) {
				print DIS "$line\n";
			}
			elsif ($chr eq '*') {
				next;
			}
			else {
				if (&check_end($chr, $final, \%Chrm) == 1) {
					print OUT "$line\n";
				}
				else {
					print DIS "$line\n";
					$bad_id{$data[0]} = 1;
				}
			}
		}
	}
}
close OUT;
close DIS;

sub get_length($$$) {
	my ($pos, $length, $cygar) = @_;
	my @info = split (//, $cygar);
	my ($tot_M, $tot_H, $tot_S, $tot_D, $tot_I, $n) = (0, 0, 0, 0, 0, 0);
	my @number;
	foreach my $digit (@info) {
		#print "$digit\t";
		if ($digit =~ m/\d+/) {
			$number[$n] = $digit;
			$n++;
		}
		elsif ($digit eq 'M') {
			my $value = join ('', @number);
			$tot_M += $value;
			$n = 0;
			undef @number;
		}
		elsif ($digit eq 'H') {
			my $value = join ('', @number);
			$tot_H += $value;
			$n = 0;
			undef @number;
		}
		elsif ($digit eq 'S') {
			my $value = join ('', @number);
			$tot_S += $value;
			$n = 0;
			undef @number;
		}
		elsif ($digit eq 'D') {
			my $value = join ('', @number);
			$tot_D += $value;
			$n = 0;
			undef @number;
		}
		elsif ($digit eq 'I') {
			my $value = join ('', @number);
			$tot_I += $value;
			$n = 0;
			undef @number;
		}
	}
	#print "\n$cygar\t$tot_M + $tot_H + $tot_S + $tot_D - $tot_I\n";
	my $cygar_length = $tot_M + $tot_H + $tot_S + $tot_D - $tot_I;
	if ($length >= $cygar_length) {
		my $final = $pos + $length;
		return $final;
	}
	else {
		my $final = $pos + $cygar_length;
		return $final;
	}
}

sub get_chr ($$$) {
	my ($line, $chr, $end) = @_;
	my @data = split (/\s+/, $line);
	foreach my $field (@data) {
		if ($field =~ m/^SN:/) { $$chr = (split(/:/, $field))[1]; }
		elsif ($field =~ m/^LN:/) {$$end = (split(/:/, $field))[1];}
	}
}

sub check_end ($$$) {
	my ($chr, $end, $Chrm) = @_;
	if ($end > $$Chrm{$chr}{end}) {
		return 0;
	}
	else {
		return 1;
	}
}















