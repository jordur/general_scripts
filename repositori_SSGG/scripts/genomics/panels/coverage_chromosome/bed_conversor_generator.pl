#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (sum);
use File::Basename;
use Getopt::Long;
use Getopt::Std;

my $input_csv_gene = $ARGV[0];
my $input_gene_directory= $ARGV[1];
my $output_bed_gene_directory= $ARGV[2];

my @gene_list = &get_gene_list($input_csv_gene);
## Get gene list


&convert_bed($input_gene_directory, $output_bed_gene_directory, @gene_list);
`chmod +x $output_bed_gene_directory/convert_commands.sh`;
system("$output_bed_gene_directory/convert_commands.sh");
`rm $output_bed_gene_directory/convert_commands.sh`;

sub convert_bed{
	my ($input_gene_directory,$output_bed_gene_directory, @gene_list) = @_;
	open CONVERT_COMMANDS, ">$output_bed_gene_directory/convert_commands.sh";
	
	foreach my $gene (@gene_list){
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_totals_coords.txt -o $output_bed_gene_directory/$gene"."_totals_coords.bed -ti chr-coord -to bed\n";
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_exon_coords.txt -o $output_bed_gene_directory/$gene"."_exon_coords.bed -ti chr-coord -to bed\n";
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_splicing_coords.txt -o $output_bed_gene_directory/$gene"."_splicing_coords.bed -ti chr-coord -to bed\n";
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_strong_splicing_coords.txt -o $output_bed_gene_directory/$gene"."_strong_splicing_coords.bed -ti chr-coord -to bed\n";
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_UTR_coords.txt -o $output_bed_gene_directory/$gene"."_UTR_coords.bed -ti chr-coord -to bed\n";
		print CONVERT_COMMANDS "transform_intervals_files.pl -i $input_gene_directory/$gene"."_vars_coords.txt -o $output_bed_gene_directory/$gene"."_vars_coords.bed -ti chr-coord -to bed\n";
	}
	close CONVERT_COMMANDS;
}

sub get_gene_list{
	my $input_csv_gene = $_[0];
	my @gene_list;
	
	open GENE_CSV, "$input_csv_gene" or die "ERROR: Couldn't open CSV gene file !";
	while (<GENE_CSV>){
		chomp $_;
		push (@gene_list, $_); 		
	}
	close GENE_CSV;
	
		
	return @gene_list;
}