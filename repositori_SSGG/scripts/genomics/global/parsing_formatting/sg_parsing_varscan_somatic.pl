#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic VarScan in target regions
# @Author: JM Rosa
# @Contributors:
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
use FileHandle;

sub main ($);
sub check_vars ($);
sub get_alleles ($);
sub get_info ($$);
sub get_type ($$);
sub get_SB ($$);
sub get_multiple_alleles ($);
sub get_data_multiple ($$);
sub get_info_single ($$$);
sub get_value ($);

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

sub main ($) {
	my $config = shift;
	my $log_file_handle = $config->{log_file_handle}; 
	my $out_file_handle = $config->{out_file_handle};
	my $error_file_handle = $config->{error_file_handle};
	my $input = $config->{in_file_handle};

	while (my $line = <$input>) {
		chomp $line;
		my @data = split (/\s+/, $line);

		my $check = &check_vars (\@data);
		
		if ($check eq 'regular_single') {
			my $string = "$data[0]\t$data[1]\t";
			my @alleles = &get_alleles (\@data);
			$string .= "$alleles[0]\t$alleles[1]\t";
			my @freqs = &get_info (\@data, "freq");
			my @genos = &get_info (\@data, "geno");
			my @depths = &get_info (\@data, "depth");
			my $type =  &get_type ($genos[0], $genos[1]); 
			$string .= "$type\t$depths[1]\t$genos[1]\t$freqs[1]\t$depths[0]\t$genos[0]\t$freqs[0]\t";
			my $sb = &get_SB ($data[5], $line);
			$string .= "$sb\n";
			print $out_file_handle $string;
		}
		elsif ($check eq 'regular_multiple') {
			my @alleles = &get_multiple_alleles (\@data);
			for (my $i = 0; $i < scalar (@alleles); $i++) {
				my $string = "$data[0]\t$data[1]\t";
				$string .= "$alleles[$i]{ref}\t$alleles[$i]{var}\t";
				my @info = &get_data_multiple (\@data, \%{$alleles[$i]});
				my @freqs = &get_info (\@info, "freq");
				my @genos = &get_info (\@info, "geno");
				my @depths = &get_info (\@info, "depth");
				my $type =  &get_type ($genos[0], $genos[1]); 
				$string .= "$type\t$depths[1]\t$genos[1]\t$freqs[1]\t$depths[0]\t$genos[0]\t$freqs[0]\t";
				my $sb = &get_SB ($info[5], $line);
				$string .= "$sb\n";
				print $out_file_handle $string;
			}
		}
		elsif ($check eq 'no_depth_normal_single') {
			my $string = "$data[0]\t$data[1]\t";
			my @alleles = &get_alleles (\@data);
			$string .= "$alleles[0]\t$alleles[1]\t";
			my $freqs = &get_info_single (\@data, "freq", "11");
			my $genos = &get_info_single (\@data, "geno", "11");
			my $depths = &get_info_single (\@data, "depth", "11");
			my $type = "unknown"; 
			$string .= "$type\t$depths\t$genos\t$freqs\t-\t-\t-\t";
			my $sb = &get_SB ($data[5], $line);
			$string .= "$sb\n";
			print $out_file_handle $string;
		}
		elsif ($check eq 'no_depth_tumour_single') {
			my $string = "$data[0]\t$data[1]\t";
			my @alleles = &get_alleles (\@data);
			$string .= "$alleles[0]\t$alleles[1]\t";
			my $freqs = &get_info_single (\@data, "freq", "10");
			my $genos = &get_info_single (\@data, "geno", "10");
			my $depths = &get_info_single (\@data, "depth", "10");
			my $type = "unknown"; 
			$string .= "$type\t-\t-\t-\t$depths\t$genos\t$freqs\t";
			my $sb = &get_SB ($data[5], $line);
			$string .= "$sb\n";
			print $out_file_handle $string;
		}
		else {
			print $error_file_handle "$line\n";
		}
	}

}

sub get_info_single ($$$) {
	my ($data, $info, $col) = @_;
	my @info = split (/:/, $$data[$col]);
	my $freq = sprintf('%.2f', ($info[3] / $info[1]));

	if ($info eq 'freq') {
		return $freq;
	}
	elsif ($info eq 'geno') {
		return &get_value ($freq);
	}
	elsif ($info eq 'depth') {
		return $info[1];
	}
}

sub get_data_multiple ($$) {
	my ($data, $allele) = @_;
	my @data;
	@data = @{$data};
	if ($$allele{geno} =~ m/[-\+]/) {
		my $string = '\\'.$$allele{geno}.":";
		if ($data[10] !~ m/$string/) {
			my $depth = (&get_info (\@data, "depth"))[1];
			my @info = split (/:/, $data[10]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[10] = join (':', @info);
		}
		elsif ($data[11] !~ m/$string/) {
			my $depth = (&get_info (\@data, "depth"))[1];
			my @info = split (/:/, $data[11]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[11] = join (':', @info);
		}
	}
	else {
		my $string = $$allele{var}.":";
		my $string_1 = $$allele{geno}.":";
		if ($data[10] !~ m/$string/ and $data[10] !~ m/$string_1/) {
			my $depth = (&get_info (\@data, "depth"))[1];
			my @info = split (/:/, $data[10]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[10] = join (':', @info);
		}
		elsif ($data[11] !~ m/$string/ and $data[11] !~ m/$string_1/) {
			my $depth = (&get_info (\@data, "depth"))[1];
			my @info = split (/:/, $data[11]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[11] = join (':', @info);
		}
	}
	return @data;
}

sub get_multiple_alleles ($) {
	my $data = shift;
	my @alleles;
	my $ref = $$data[2];
	my $n = 0;
	my @vars = split(/\//, $$data[3]);
	foreach my $var (@vars) {
		if ($var =~ m/-/) {
			$alleles[$n]{ref} = $ref.$var;
			$alleles[$n]{ref} =~ s/-//;
			$alleles[$n]{var} = $ref;
			$alleles[$n]{geno} = $var; 
		}		
		elsif ($var =~ m/\+/) {
			$alleles[$n]{ref} = $ref;
			$alleles[$n]{var} = $ref.$var;
			$alleles[$n]{var} =~ s/\+//;
			$alleles[$n]{geno} = $var; 
		}		
		else {
			$alleles[$n]{ref} = $ref;
			$alleles[$n]{var} = $var;
			$alleles[$n]{geno} = &get_symbol ($ref, $var); 
		}
		$n++;
	}
	#print "@alleles\n";
	return @alleles;
}

sub get_symbol ($$) {
	my ($ref, $var) = @_;
		my %symbols;
		$symbols{'A'}{'G'} = 'R';
		$symbols{'G'}{'A'} = 'R';
		$symbols{'C'}{'T'} = 'Y';
		$symbols{'T'}{'C'} = 'Y';
		$symbols{'G'}{'C'} = 'S';
		$symbols{'C'}{'G'} = 'S';
		$symbols{'A'}{'T'} = 'W';
		$symbols{'T'}{'A'} = 'W';
		$symbols{'G'}{'T'} = 'K';
		$symbols{'T'}{'G'} = 'K';
		$symbols{'A'}{'C'} = 'M';
		$symbols{'C'}{'A'} = 'M';
		
		return $symbols{$ref}{$var};
}

sub get_info ($$) {
	my ($data, $info) = @_;
	my @result;
	my @normal_info = split (/:/, $$data[10]);
	my @tumour_info = split (/:/, $$data[11]);
	my $normal_freq = sprintf('%.2f', ($normal_info[3] / $normal_info[1]));
	my $tumour_freq = sprintf('%.2f', ($tumour_info[3] / $tumour_info[1]));

	if ($info eq 'freq') {
		$result[0] = $normal_freq;
		$result[1] = $tumour_freq;
	}
	elsif ($info eq 'geno') {
		$result[0] = &get_value ($normal_freq);
		$result[1] = &get_value ($tumour_freq);
	}
	elsif ($info eq 'depth') {
		$result[0] = $normal_info[1];
		$result[1] = $tumour_info[1];
	}
	return @result;
}

sub get_alleles ($) {
	my $data = shift;
	my @alleles;
	if ($$data[3] =~ m/-/) {
		$alleles[0] = $$data[2].$$data[3];
		$alleles[0] =~ s/-//;
		$alleles[1] = $$data[2];
	}
	elsif ($$data[3] =~ m/\+/) {
		$alleles[1] = $$data[2].$$data[3];
		$alleles[1] =~ s/\+//;
		$alleles[0] = $$data[2];
	}
	else {
		$alleles[0] = $$data[2];
		$alleles[1] = $$data[3];
	}
	return @alleles;
}

sub get_SB ($$) {
	my ($data, $line) = @_;
	my @info = split (/:/, $data);
	if (scalar (@info) > 6) {
		print_out ($line, "ERR");
		return 'NA';
	}
	else {
		my $sb;
		my $R1_minus = $info[1];
		my $R1_plus = $info[2];
		my $R2_minus = $info[3];
		my $R2_plus = $info[4];
		
		my $minus_reads = $R1_minus + $R2_minus;
		my $plus_reads = $R1_plus + $R2_plus;
		
		if ($minus_reads == 0) { 
			$sb = '+++';
		}
		elsif ($plus_reads == 0) {
			$sb = '---';
		}
		else {
			$sb = sprintf('%.2f', ($plus_reads/$minus_reads));
		}
		return $sb;
	}
}

sub get_type ($$) {
	my ($normal_gen, $tumour_gen) = @_;
	if ($normal_gen eq $tumour_gen) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/ref/i and $tumour_gen =~ m/ref/i) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/ref/i and $tumour_gen !~ m/ref/i) {
		return "Somatic";
	}
	elsif ($normal_gen !~ m/ref/i and $tumour_gen =~ m/ref/i) {
		return "LOH";
	}
	elsif ($normal_gen !~ m/var/i and $tumour_gen =~ m/var/i) {
		return "LOH";
	}
	elsif ($normal_gen =~ m/hetero/i and $tumour_gen =~ m/hetero/i) {
		return "Germline";
	}
	elsif ($normal_gen =~ m/var/i and $tumour_gen =~ m/hetero/i) {
		return "GOH";
	}
	else {
		return "unknown";
	}
}

sub get_value ($) {
	my $freq = shift;
	my $ratio = sprintf('%.2f', $freq);
	my $value;
	if ($ratio == 0) {
		$value = "Homo_ref";
	}
	elsif ($ratio > 0 && $ratio <= 0.12) {
		$value = "P_Homo_ref";
	}
	elsif ($ratio > 0.12 && $ratio < 0.35) {
		$value = "UNC_Hetero";
	}
	elsif ($ratio >= 0.35 && $ratio < 0.65) {
		$value = "P_Hetero";
	}
	elsif ($ratio >= 0.65 && $ratio < 0.85) {
		$value = "UNC_Homo";
	}
	else {
		$value = "P_Homo_var";
	}
	return $value;
}

sub check_vars ($) {
	my $data = shift;
	my $check = "NA";
	
	if ($$data[10] !~ m/^N/ and $$data[11] !~ m/^N/) {
		$check = "regular_";
	}
	elsif ($$data[10] =~ m/^N/ and $$data[11] !~ m/^N/) {
		$check = "no_depth_normal_";
	}
	elsif ($$data[10] !~ m/^N/ and $$data[11] =~ m/^N/) {
		$check = "no_depth_tumour_";
	}
	
	if ($$data[3] !~ m/\//) {
		$check .= "single";
	}
	elsif ($$data[3] =~ m/\//) {
		$check .= "multiple";
	}
	
	return $check;
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}

sub get_log_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $log_file_handle = new FileHandle;
	$log_file_handle->open(">".$config->{log_file}) or die("ERROR: Could not write to output file ", $config->{log_file}, "\n");
	
	# make header
	my $path = `pwd`;
	chomp $path;
	print $log_file_handle "sg_parsing_varscan_somatic.pl -i $config->{input_file} -n $config->{normal} -t $config->{tumour} -o $config->{output_file} in $path\n";
	print $log_file_handle "\nStarting process... ".localtime()."\n";
	
	return $log_file_handle;
}

sub get_error_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $error_file_handle = new FileHandle;
	$error_file_handle->open(">".$config->{error_file}) or die("ERROR: Could not write to output file ", $config->{error_file}, "\n");
	
	# make header

	return $error_file_handle;
}

sub get_out_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$config->{output_id}) or die("ERROR: Could not write to output file ", $config->{output_id}, "\n");
	
	# make header
	print $out_file_handle "#Chr\tPosition\tRef_Allele\tVar_Allele\tType\t$config->{tumour}\_depth\t$config->{tumour}\_geno\t$config->{tumour}\_freq\t$config->{normal}\_depth\t$config->{normal}\_geno\t$config->{normal}\_freq\tSB\n";
	
	return $out_file_handle;
}

sub configure {
	my $args = shift;
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		
		# name options
		'normal=s',
		'tumour=s',
		
		# output options
		'output_id=s',
		'verbose',
		'quiet',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#------------------------------------------------#
# Parsing Somatic Vars from VarScan script v 1.0 #
#------------------------------------------------#

version 1.0

By Juan Manuel Rosa @ Sistemas Genomicos

Configuration options:

INTRO
		print $header;
	}	
	# set defaults
	
	my $input = $config->{input_file};
	$input =~ s/\.\.\/+//g;
	print "$config->{input_file} -> $input\n";
	my @output_name = split (/\./, $input); 
	print "$config->{input_file} -> $output_name[0]\n";
	$config->{output_file} ||= $output_name[0]."_somatic";
	$config->{output_id}  = $config->{output_file}.".snvs";
	$config->{log_file} = $config->{output_file}.".log";
	$config->{error_file} = $config->{output_file}.".error";
	
	# get input file handles
	$config->{in_file_handle} = &get_in_file_handle($config);
	
	# get error file handles
	$config->{error_file_handle} = &get_error_file_handle($config);	
	
	# get log file handles
	$config->{log_file_handle} = &get_log_file_handle($config);	
	
	# configure output file
	$config->{out_file_handle} = &get_out_file_handle($config);
	
	return $config;
}

sub usage {
	my $usage =<<END;
#------------------------------------------------#
# Parsing Somatic Vars from VarScan script v 1.0 #
#------------------------------------------------#

version 1.0

By Juan Manuel Rosa @ Sistemas Genomicos

Usage:
sg_parsing_varscan_somatic.pl [arguments]

Options
=======

--help	 			Display this message and quit

-i | --input_file	Input file from VarScan mpileup2cns. Normal sample must be first

-o | --output_id	Output name. Optional 
                       
--normal	Normal 	Sample Name
--tumour	Tumour 	Sample Name

END

	print $usage."\n".localtime()."\n";
}

exit;
