#!/usr/bin/env perl

=head1 LICENSE

  This script was developed by Sistemas Genomicos SL

=head1 CONTACT

  Please email comments or questions to:
  jm.rosa@sistemasgenomicos.com

=cut

=head1 NAME

Parsing Freqs and Populations from VEP

Version 1.0

by Juan Manuel Rosa (jm.rosa@sistemasgenomicos.com)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Parallel::ForkManager;
use Data::Dumper;


sub debug ($);
sub get_out_file_handle ($);
sub getting_annotation ($$);
sub getting_completation ($$);


# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine

&main($config);

# this is the main sub-routine - it needs the configured $config hash

#### OJO con el numero de muestras, hay que parsear columnas!!!!

sub main {
	my $n = 0;
	my $config = shift;
	my @chrs = split (/,/, $config->{chr});
	my %pops;
	my $output_file = &get_out_file_handle("$config->{out_file}.csv", $chrs[0]);
	print $output_file join '|', &get_header($chrs[0]);
	print $output_file "\n";
	my $n_samples = &get_samples($chrs[0]);
	my $rs_col = 5+($n_samples*3)+3;
	my $freq_col = 5+($n_samples*3)+4;
	my $pop_col = 5+($n_samples*3)+18;
	my $pop_file = &get_out_file_handle("$config->{out_file}.pop", $chrs[0]);
	my $chrs_count = 0;
	foreach my $chr (@chrs) {
		# get input file handles
		my $input = &get_in_file_handle($config, $chr);
		if (defined($input)){
			$chrs_count ++;
			while (my $line = <$input>){
				chomp $line;
				my @data = split (/\s+/, $line);
				if ($data[0] =~ m/^#HGNC_ID/) {
					next;
				}	
				elsif ($data[0] =~ m/^#Code/ or $data[0] =~ m/^\d+$/ or $data[0] eq '' ){
					next;
				}
				# Filtering in this script is about to be deprecated!!
				else {
#					if (defined ($config->{filtering})) {
#						my $cons_col = 5+($n_samples*3)+7;
#						if ($data[$cons_col] =~ m/NMD_TRANS/g) {next;}
#						if ($data[$cons_col] =~ m/WITHIN_NON/g) {next;}
#						if ($data[$cons_col] eq 'INTRONIC') {next;}
#						if ($data[$cons_col] eq 'DOWNSTREAM') {next;}
#						if ($data[$cons_col] eq 'UPSTREAM') {next;}
#						
#						if ($data[$freq_col] eq '-') {
#							print $output_file join ('|', @data),"\n";
#						}
#						else {
#							$data[$freq_col] = &parse_freq ($data[$freq_col]);
#							if (!defined $pops{vars}{$data[$rs_col]}) {
#								my $pops = &check_pops($config, $data[$pop_col], \%pops, \$n, $chr);
#								$pops{vars}{$data[$rs_col]}{pops} = $pops;
#							}
#							$data[$pop_col] = $pops{vars}{$data[$rs_col]}{pops};
#							print $output_file join ('|', @data),"\n";
#						}	
#					}
#					else {
						if ($data[$freq_col] eq '-') {
							print $output_file join ('|', @data),"\n";
						}
						else {
							$data[$freq_col] = &parse_freq ($data[$freq_col]);
							if (!defined $pops{vars}{$data[$rs_col]}) {
								my $pops = &check_pops($config, $data[$pop_col], \%pops, \$n, $chr);
								$pops{vars}{$data[$rs_col]}{pops} = $pops;
							}
							$data[$pop_col] = $pops{vars}{$data[$rs_col]}{pops};
							print $output_file join ('|', @data),"\n";
						}	
#					}	
				}
			}	
		}
	}
	if ($chrs_count == 0){
		die "ERROR: Either no variants were found for any of the chromosomes, or no input files were found!!";
	}
	print $pop_file "#Code\tPop_Name\n";
	for (my $i=1; $i <= $n; $i++){
		print $pop_file "$i\t$pops{populations}{$i}{name}\n";
	}
}

debug("Finished!") if defined $config->{verbose};
exit;

sub parse_freq {
	my $freqs = shift;
	my @freqs = split (/\//, $freqs);
	my $string = '';
	foreach my $freq (@freqs) {
		$string .= '['.$freq.']';
	}
	return $string;
}

sub check_pops {
	my ($config, $data, $pops, $n, $chr) = @_;
	my @alleles = split (/\//, $data);
	my $string;
	my $file =  $config->{input_file}.".$chr.snvs";
	foreach my $allele (@alleles) {
		$string .= '[';
		my @pops = split (/,/, $allele);
		foreach my $pop (@pops){
			# If no information for populations is available in the annotation file:
			if ($pop eq '-' or $pop eq ''){
				print "INFO: Non defined allele $allele from chr$chr and data $data\n";
			} else {
				# If information is available:
				my $name = `awk \'\$1==$pop \{print \$2\}\' $file`;
				chomp $name;
				if (defined $$pops{populations}{$name}){
					$string .= "$$pops{populations}{$name}{number},";
				} else {
					$$n++;
					$$pops{populations}{$name}{number} = $$n;
					$string .= "$$pops{populations}{$name}{number},";
					$$pops{populations}{$$n}{name} = $name;
				}
			}
		}
		$string =~ s/,$//;
		$string .= "\]";	
	}
	return $string;
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
		
		# global options
		'pops=s',
		'chr=s',
		
		# output options
		'output_file=s',
		'filtering',
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
#---------------------------------------------#
# PARSING FREQS AND POPULATIONS FROM ENSEMBL  #
#---------------------------------------------#

version 1.0

By sistemas Genomicos

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	# set defaults
	
	$config->{pop} ||= "all";
	$config->{chr} ||= "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y";
	
	my $input = $config->{input_file};
	$input =~ s/\.\.\/+//g;
	my @output_name = split (/\./, $input); 
	$config->{output_file} ||= $output_name[0]."_parsed_pops";
	$config->{out_file} = $config->{output_file};
	$config->{log_file} = $config->{output_file}.".log";

	return $config;
}

# gets file handle for input
sub get_in_file_handle {
	my $config = shift;
	my $chr = shift;

	# define the filehandle to read input from
	my $in_file_handle;
	
	if(defined($config->{input_file})) {    
		# check defined input file exists
		unless (-e "$config->{input_file}.$chr.snvs"){
			print("WARNING: Could not find input file ", $config->{input_file}.".$chr.snvs. Please check that corresponding chromosome is not included in design!!", "\n")
		} else{
			$in_file_handle = new FileHandle;
			$in_file_handle->open( "$config->{input_file}.$chr.snvs" ) or die("ERROR: Could not read from input file $config->{input_file}.$chr.snvs", "\n");
        }
    }
	# no file specified - try to read data off command line
	else {
		$in_file_handle = 'STDIN';
		debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
	}
	return $in_file_handle;
}

#get file for output
sub get_out_file_handle ($) {
	my $file = shift;
	my $chr = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$file) or die("ERROR: Could not write to output file ", $file, "\n");

	return $out_file_handle;
}

sub get_header {
	my $chr = shift;
	my $line = `grep -m 1 "^#HG" $config->{input_file}.$chr.snvs`;
	chomp $line;
	my @results = split (/\s+/, $line);
	return @results;	
}

sub get_samples {
	my $chr = shift;
	my $line = `grep -m 1 "^#HG" $config->{input_file}.$chr.snvs`;
	chomp $line;
	my $number =()= $line =~ /ratio/gi;
	#print "$number\n";
	return $number;	
}

sub usage {
	my $usage =<<END;
#---------------------------------------------#
# PARSING FREQS AND POPULATIONS FROM ENSEMBL  #
#---------------------------------------------#

version 1.0

By Sistemas Genomicos

Usage:
perl sg_parsing_freqs_pops.pl [arguments]

Options
=======

--help                 Display this message and quit

-i | --input_file      Input files root name

-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet 

--filtering		To filter by consequences (NMD_TRANS, WITHIN_NON, INTRONIC, DOWNSTREAM, UPSTREAM)					 
--pop		   Population names (comma separated) to filter

--chr			Chromosome list (comma separated, default 'all no MT')

END

	print $usage;
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}

# prints debug output with time
sub debug ($) {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}