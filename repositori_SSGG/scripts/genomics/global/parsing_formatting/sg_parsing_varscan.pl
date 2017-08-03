#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing VarScan in target regions
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
use Math::BigFloat;
use FileHandle;
use genomics::genotype;

sub main ($);
sub check_vars ($$);
sub get_alleles ($);
sub get_info ($$$$);
sub get_SB ($$$);
sub get_multiple_alleles ($);
sub get_data_multiple ($$$);
sub get_info_single ($$$);

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
		if ($data[0] eq 'Chrom') { next;}
		
		my $check = &check_vars ($config, \@data);
		
		if ($check eq 'na') {
			next;
		}		
		elsif ($check eq 'single') {
			my $string = "$data[0]\t$data[1]\t";
			my @alleles = &get_alleles (\@data); ############ modificar para coger indels
			$string .= "$alleles[0]\t$alleles[1]\t";
			my $n = 10;
			my @samples = split(/,/, $config->{sample});
			my $type = '-'; 
			$string .= "$type\t";
			foreach my $sample (@samples) {
				my ($freq, $geno, $depth);
				if ($data[$n] =~ m/^N/) {
					$freq = '-';
					$geno = '-';
					$depth = '-';
				}
				else {
					$freq = &get_info (\@data, "freq", $n, 2);
					$geno = &get_info (\@data, "geno", $n, 2);
					$depth = &get_info (\@data, "depth", $n, 2);
				}
				$string .= "$depth\t$geno\t$freq\t";
				$n++;
			}
			my $sb = &get_SB ($data[5], $line, $error_file_handle);
			$string .= "$sb\n";
			print $out_file_handle $string;
		}
		elsif ($check eq 'multiple') {
			my @indel;
			my @alleles = &get_multiple_alleles (\@data, \@indel); ##################
			for (my $i = 0; $i < scalar (@alleles); $i++) {
				my $string = "$data[0]\t$data[1]\t";
				$string .= "$alleles[$i]{ref}\t$alleles[$i]{var}\t";
				my $n = 10;
				my @samples = split(/,/, $config->{sample});
				my $type = '-'; 
				$string .= "$type\t";
				foreach my $sample (@samples) {
					my ($freq, $geno, $depth);
					if ($data[$n] =~ m/^N/) {
						$freq = '-';
						$geno = '-';
						$depth = '-';
					}
					else {
						my @info = &get_data_multiple (\@data, \%{$alleles[$i]}, $n);
						$freq = &get_info (\@info, "freq", $n, 0);
						$geno = &get_info (\@info, "geno", $n, $indel[$i]);
						$depth = &get_info (\@info, "depth", $n, 0);
					}
					$string .= "$depth\t$geno\t$freq\t";
					$n++;
				}
				my $sb = &get_SB ($data[5], $line, $error_file_handle);
				$string .= "$sb\n";
				print $out_file_handle $string;
			}
		}
		else {
			print $error_file_handle "$line\n";
		}
	}

}

sub get_data_multiple ($$$) {
	my ($data, $allele, $pos) = @_;
	my @data;
	@data = @{$data};
	if ($$allele{geno} =~ m/\+/ or $$allele{geno} =~ m/-/) {
		my $string = '\\'.$$allele{geno}.":";
		if ($data[$pos] !~ m/$string/) {
			my $depth = &get_info (\@data, "depth", $pos, 0);
			my @info = split (/:/, $data[$pos]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[$pos] = join (':', @info);
		}
	}
	else {
		my $string = $$allele{var}.":";
		my $string_1 = $$allele{geno}.":";
		if ($data[$pos] !~ m/$string/ and $data[$pos] !~ m/$string_1/) {
			my $depth = &get_info (\@data, "depth", $pos, 0);
			my @info = split (/:/, $data[$pos]);
			$info[2] = $info[1] - $depth;
			$info[3] = '0';
			$data[$pos] = join (':', @info);
		}
	}
	return @data;
}

sub get_multiple_alleles ($) {
	my ($data, $indel) = @_;
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
			$$indel[$n] = 1;
		}
		elsif ($var =~ m/\+/) {
			$alleles[$n]{ref} = $ref;
			$alleles[$n]{var} = $ref.$var;
			$alleles[$n]{var} =~ s/\+//;
			$alleles[$n]{geno} = $var;
			$$indel[$n] = 1;
		}
		else {
			$alleles[$n]{ref} = $ref;
			$alleles[$n]{var} = $var;
			$alleles[$n]{geno} = &get_symbol ($ref, $var); 
			$$indel[$n] = 0;
		}
		$n++;
	}
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

sub get_info ($$$$) {
	my ($data, $info, $pos, $indel) = @_;
	my $result;
	my $ref = $$data[2];
	my $var;
	if ($indel == 1){
		# This creates an fictitious indel
		$var = $ref . $ref;
	} elsif ($indel == 0){
		# This assures that the variant isn't an indel
		$var = $ref;
	} else {
		# If only one variant allele available
		$var = $$data[3];
	}
	my @info = split (/:/, $$data[$pos]);
	my $freq = sprintf('%.2f', (substr($info[4],0,length($info[4])-1) / 100)); # sprintf('%.2f', ($info[3] / $info[1]));
	
	if ($info eq 'freq') {
		$result = $freq;
	}
	elsif ($info eq 'geno') {
		my @tmp = split(",", genotype->get_genotype ($info[1],$freq,$ref,$var));
		$result = $tmp[0];
	}
	elsif ($info eq 'depth') {
		$result = $info[1];
	}
	return $result;
}

sub get_alleles ($) {
	my ($data) = @_;
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

sub get_SB ($$$) {
	my ($data, $line, $error_file_handle) = @_;
	my @info = split (/:/, $data);
	if (scalar (@info) > 6) {
		print $error_file_handle "$line\n"; 
		return 'NA';
	}
	else {
		# Ref are alleles for the reference, var for the variant, and minus stays for negative strand, while plus for positive strand
		my $ref_minus = $info[1];
		my $ref_plus = $info[2];
		my $ref_all = Math::BigFloat->new($ref_minus + $ref_plus);
		my $var_minus = $info[3];
		my $var_plus = $info[4];
		my $var_all = Math::BigFloat->new($var_minus + $var_plus);
		my $plus_all = $ref_plus + $var_plus;
		my $all = Math::BigFloat->new($ref_all + $var_all);
		
		# Phred scale p-value from Fisher's exact test is applied to reference and variant alleles count over positive and negative strands:
		my $sb = sprintf "%.3f", -10*log(($ref_all->bnok($ref_plus) * $var_all->bnok($var_plus)) / ($all->bnok($plus_all)));
		
		# TODO: remove this part of code after some months of test of substitution code (20130610)
		#my $minus_reads = $R1_minus + $R2_minus;
		#my $plus_reads = $R1_plus + $R2_plus;
		
		#if ($minus_reads == 0) { 
		#	$sb = '+++';
		#}
		#elsif ($plus_reads == 0) {
		#	$sb = '---';
		#}
		#else {
		#	$sb = sprintf('%.2f', ($plus_reads/$minus_reads));
		#}
		return $sb;
	}
}

sub check_vars ($$) {
	my ($config, $data) = @_;
	my $check = "NA";
	my @samples = split (/,/, $config->{sample});
	my $n = 10;
	my $value = 0;
	foreach my $sample (@samples) {
		if ($$data[$n] =~ m/^N/) {
			$value++;
		}
		$n++;
	}
	if ($value == scalar(@samples)) {
		$check = 'na';
	}
	else {
		if ($$data[3] !~ m/\//) {
			$check = "single";
		}
		elsif ($$data[3] =~ m/\//) {
			$check = "multiple";
		}
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
	print $log_file_handle "sg_parsing_varscan.pl -i $config->{input_file} --sample $config->{sample} -o $config->{output_file} in $path\n";
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
	
	print $out_file_handle "#Chr\tPosition\tRef_Allele\tVar_Allele\tType\t";
	my @samples = split(/,/, $config->{sample});
	foreach my $sample (@samples) {
		print $out_file_handle "$sample\_depth\t$sample\_geno\t$sample\_freq\t";
	}
	print $out_file_handle "SB\n";
	
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
		'sample=s',
		
		# output options
		'output_file=s',
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
	#print "$config->{input_file} -> $input\n";
	my @output_name = split (/\./, $input); 
	#print "$config->{input_file} -> $output_name[0]\n";
	$config->{output_file} ||= $output_name[0]."_varscan_parsed";
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
#--------------------------------------------#
# Parsing Variants from VarScan script v 1.0 #
#--------------------------------------------#

version 1.0

By Juan Manuel Rosa @ Sistemas Genomicos

Usage:
sg_parsing_varscan.pl [arguments]

Options
=======

--help	 			Display this message and quit

-i | --input_file	Input file from VarScan mpileup2cns. 

-o | --output_file	Output name. Optional 
                       
--sample	Sample Name (comma separated for multiple)

END

	print $usage."\n".localtime()."\n";
}

exit;
