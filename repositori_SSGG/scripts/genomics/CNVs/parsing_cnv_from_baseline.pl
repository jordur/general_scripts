#!/usr/bin/env perl
###########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to parse CNV regions from baseline analysis
# @Author: JM Rosa
# @Contributors: Arbol
############################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use analysis_region;
use file_handle;


# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $tumour_id;
my $normal_id;
my $output_id = cwd() . "/output_file";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"t=s"	=> \$tumour_id,
				"n=s"	=> \$normal_id,
				"o=s"	=> \$output_id,
);

print "\n=================================================================
$basename: Script to parse CNV regions from baseline analysis v1.0\n";

if (not defined $tumour_id or not defined $normal_id) {
	die "
			Options:
		-t Tumour CNATable File
		-n Normal CNATable File
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	my $path = `pwd`;
	chomp $path;
	
	print "Running $basename with parameters t=$tumour_id, n=$normal_id, o=$output_id\n";
	
	print "Printing data to tmp file... ".localtime()."\n";
	
	&get_tmp_file ("$output_id", "tmp");
	
	print "Parsing data to print results... ".localtime()."\n";
	
	&parsing_results ("$output_id.tmp", "$output_id", "txt");
}

sub parsing_results {
	my ($input, $output, $extension) = @_;
	`sort -k1n -k2n $input | sed -e 's/^23/X/' -e 's/^24/Y/' -e 's/^25/M/' > $output.$extension`;
	return 1;	
}

sub get_tmp_file {
	my ($file, $extension) = @_;
	### Handling files ###	
	my $tumour_fh = &get_in_file_handle($tumour_id);
	my $normal_fh = &get_in_file_handle($normal_id);
	my $ouput = &get_out_file_handle($file, $extension);
	print $ouput &header(),"\n";
	
	### Setting variables ###
	my $compare = 'N/A';
	my $tumour_line;
	my $normal_line;
	
	### Computing values ###
	do {
		if ($compare eq 'N/A') {
			$tumour_line = <$tumour_fh>;
			$normal_line = <$normal_fh>;
			
			if ($tumour_line =~ m/^Targeted/) {
				$tumour_line = <$tumour_fh>;
			} 
			
			if ($normal_line =~ m/^Targeted/) {
				$normal_line = <$normal_fh>;
			} 
			$compare = &compare_lines ($tumour_line, $normal_line);			
		}
		elsif ($compare == 0) {
			print $ouput &print_output($compare, $tumour_line, $normal_line),"\n";

			$tumour_line = <$tumour_fh>;
			$normal_line = <$normal_fh>;

			$compare = &compare_lines ($tumour_line, $normal_line);
		}
		elsif ($compare == 1){
			print $ouput &print_output($compare, $tumour_line, $normal_line),"\n";

			$tumour_line = <$tumour_fh>;

			$compare = &compare_lines ($tumour_line, $normal_line);				
		}
		elsif ($compare == 2){
			print $ouput &print_output($compare, $tumour_line, $normal_line),"\n";

			$normal_line = <$tumour_fh>;

			$compare = &compare_lines ($tumour_line, $normal_line);						
		}
		
		else {
			die "Strange value for comparison: $compare\n";
		}		
	
	} while (!eof($tumour_fh) and !eof($normal_fh)); 
	
	return 1;
}

sub print_output {
	my ($compare, $tumour, $normal) = @_;
	my $result;
	my $t_object = analysis_region->get_region($tumour);
	my $n_object = analysis_region->get_region($normal);
	
	if ($compare == 0){
		my $value = analysis_region->get_log_ratio ($tumour, $normal);
		$result = "$t_object->{c_chr}\t$t_object->{start}\t$t_object->{end}\tchr$t_object->{chrom}:$t_object->{start}-$t_object->{end}\t$value";
	}
	elsif ($compare == 1){
		my $value = (split (/\s+/, $tumour))[7];
		$result = "$t_object->{c_chr}\t$t_object->{start}\t$t_object->{end}\tchr$t_object->{chrom}:$t_object->{start}-$t_object->{end}\t$value";
	}
	
	elsif ($compare == 2){
		my $value = (split (/\s+/, $normal))[7];
		$result = "$n_object->{c_chr}\t$n_object->{start}\t$n_object->{end}\tchr$n_object->{chrom}:$n_object->{start}-$n_object->{end}\t$value";
	}
	
	return $result;	
}

sub compare_lines {
	my ($tumour, $normal) = @_;
	
	my $compare;
	my $t_region = analysis_region->get_region($tumour);
	my $n_region = analysis_region->get_region($normal);
	
	if ($t_region->{c_chr} == $n_region->{c_chr}){
		if ($t_region->{start} == $n_region->{start}){
			$compare = 0;
		}
		elsif ($t_region->{start} > $n_region->{start}){
			$compare = 2;
		}
		else {
			$compare = 1;
		}
	}
	
	elsif ($t_region->{c_chr} > $n_region->{c_chr}){
		$compare = 2;
	}
	
	elsif ($t_region->{c_chr} < $n_region->{c_chr}){
		$compare = 1;
	}

	else {
		$compare = 4;
	}	
	return $compare;
}

sub header {
	my $result = "#Chr\tOriStCoordinate\tOriEndCoordinate\tChr:OriStCoordinate-OriEndCoordinate\tLogRatio";
	return $result;
}

exit;