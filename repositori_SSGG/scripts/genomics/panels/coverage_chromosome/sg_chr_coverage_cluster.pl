#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (sum);
use File::Basename;
use Getopt::Long;
use Getopt::Std;

### Input Data ###
#CSV with gene list
#my $input_csv_gene="/home/gmarco/Onco2/oncogeneprofile.csv";

#Gene BED file directory

#my $input_bed_not_covered = "/home/gmarco/Onco2/Onco2-NA12144_b_S7_not_covered.bed";
#my $input_gene_directory= "/home/gmarco/Onco2/ONCOGENEPROFILE_min0/genes";
#my $output_bed_gene_directory= "/home/gmarco/Onco2/ONCOGENEPROFILE_min0/genes/bed";
#my $input_bed_gene_directory = $output_bed_gene_directory;

## Get gene list
#&convert_bed($input_gene_directory, $output_bed_gene_directory, @gene_list);

main();

sub main {
	my $gene;
	my ($chr,$input_totals_coords,$input_exon_coords,$input_splicing_coords,$input_strong_splicing_coords,$input_UTR_coords,$input_vars_coords);
	
	my %options=();
	
	getopts("p:c:n:b:o:", \%options);
	
	print "\n=================================================================
	SG NOT COVERED PERCENTAGE PER CHROMOSOME\n";
	#print "\t",`date`;
	if (not defined $options{p} or not defined $options{c} or not defined $options{n} or not defined $options{b} or not defined $options{o}){
		print "\t!!! YOU MUST SPECIFY ALL INPUT FLAGS !!!\n";
		print "!!!! ALL DIRECTORIES INPUT PATHS END WITHOUT SLASH !!!!!\n";
		print "-p Panel name\n";
		print "-c CSV File containing the list of genes in the panel.\n";
		print "-n BED File of the non covered regions !\n";
		print "-b NO END SLASH BED folder where all the chromosome BED files of the panel are located!! ie /home/panel/genes\n";
		print "-o Output path USING ABSOLUTE PATH !!! ie:/home/resulsts/output/\n";
		print "=================================================================\n\n";
		exit;
	}
	print "=================================================================\n";

	#Check input csv file exists
	if (not -e $options{c}){
		die "ERROR: Couldn't read CSV gene list file.\n";
	}
	
	#Check input not covered BED file exists
	if (not -e $options{n}){
		die "ERROR: Couldn't read not covered regions BED file.\n";
	}
	
	#Check if output chromosome bed folder exists
	my $input_bed_gene_directory = $options{b};
	if (not defined $options{b}){
		die "ERROR: You must specify the chromosome BED folder. ie: -bf /home/panel/bed\n\n";
	}
	if (not -d $input_bed_gene_directory){
		die "ERROR: Chromosome BED input doesn't exist !!! ie: -bf /home/panel/bed\n\n";
	}
	
	#Check an output path is given
	if (not defined $options{o}){
		die "ERROR: You must specify the full output: -o /home/user/output/\n\n";
	}

	#Check if output directory exists
	my $output_dir = $options{o};
	if (not -d $output_dir){
		print "WARNING: Couldn't find output directory !\n";
		print "Creating output directory...\n";
		system ("mkdir -p $output_dir");
		if (not -d $output_dir){
			die "ERROR: Couldn't create output directory. Path incorrect or insufficient permissions.\n";
		}
	}
	
	#Input handlers
	my $panel_name = $options{p};
	my $input_csv_gene = $options{c};
	my $input_bed_not_covered = $options{n};
	
	my @gene_list = &get_gene_list($input_csv_gene);
	
	## Open output_file
	open OUTPUT, ">$output_dir/$panel_name"."_chr_not_covered_stats";
	print OUTPUT "HGNC\tCHR\tNC_TOTAL_GENE\tNC_EXON\tNC_SPLICING\tNC_STRONG_SPLICING\tNC_UTR\n";
	
	foreach $gene (@gene_list){
		
		chomp $gene;
 		
		$input_totals_coords = $input_bed_gene_directory."/$gene"."_totals_coords.bed";
		$input_exon_coords = $input_bed_gene_directory."/$gene"."_exon_coords.bed";
		$input_splicing_coords = $input_bed_gene_directory."/$gene"."_splicing_coords.bed";
		$input_strong_splicing_coords = $input_bed_gene_directory."/$gene"."_strong_splicing_coords.bed";
		$input_UTR_coords = $input_bed_gene_directory."/$gene"."_UTR_coords.bed";
#		$input_vars_coords = $input_bed_gene_directory."/$gene"."_vars_coords.bed";
		
		$chr = (split ("\n",`awk '{print \$1}' $input_totals_coords`))[0];
		
		#my $substract_command = "subtractBed -a $input_totals_coords -b ";
		my $awk_command = " | awk '{print \$3-\$2+1}'";
		my $awk = "awk '{print \$3-\$2+1}'";
		
		my $sum_total_gen = &col_sum(`$awk $input_totals_coords`);
		my $sum_total_exon_coords = &col_sum(`$awk $input_exon_coords`);
		my $sum_total_splicing_coords = &col_sum(`$awk $input_splicing_coords`);
		my $sum_total_strong_splicing_coords = &col_sum(`$awk $input_strong_splicing_coords`);
		my $sum_total_UTR_coords = &col_sum(`$awk $input_UTR_coords`);
		
		
		my $sub_total_gen = &col_sum(`subtractBed -a $input_totals_coords -b $input_bed_not_covered $awk_command`);
		my $sub_total_exon_coords = &col_sum(`subtractBed -a $input_exon_coords -b $input_bed_not_covered $awk_command`);
		my $sub_total_splicing_coords = &col_sum(`subtractBed -a $input_splicing_coords -b $input_bed_not_covered $awk_command`);
		my $sub_total_strong_splicing_coords = &col_sum(`subtractBed -a $input_strong_splicing_coords -b $input_bed_not_covered $awk_command`); 
		my $sub_total_UTR_coords = &col_sum(`subtractBed -a $input_UTR_coords -b $input_bed_not_covered $awk_command`);
#		my $sub_perc_vars;
#		
#		if (-z "$input_vars_coords") {
#			my $total_vars_coords = 0;
#			$sub_perc_vars = 0;
#		}
#		else {
#			my $total_vars_coords =  &col_sum(`subtractBed -a $input_vars_coords -b $input_bed_not_covered $awk_command`);
#			$perc_vars = &percentage_coverage($sub_total_gen, $total_vars_coords);
#		}

		my $perc_gen = &percentage_coverage($sum_total_gen, $sub_total_gen);
		my $perc_exon = &percentage_coverage($sum_total_exon_coords, $sub_total_exon_coords);
		my $perc_splicing = &percentage_coverage($sum_total_splicing_coords, $sub_total_splicing_coords);
		my $perc_strong_splicing = &percentage_coverage($sum_total_strong_splicing_coords, $sub_total_strong_splicing_coords);
		my $perc_UTR = &percentage_coverage($sum_total_UTR_coords, $sub_total_UTR_coords);
		
		print OUTPUT "$gene\t$chr\t$perc_gen\t$perc_exon\t$perc_splicing\t$perc_strong_splicing\t$perc_UTR\n";
	}
	close OUTPUT;
	
	print "DONE !! Output result file: $output_dir/$panel_name"."_chr_not_covered_stats\n\n";
	 
}

sub percentage_coverage {
	my ($total_elemento, $total_subtract) = @_;
	my $resta = $total_elemento - $total_subtract;
	my $porc = ($resta/$total_elemento)*100;
	$porc = sprintf("%.2f", $porc);
	return $porc;
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

sub col_sum {
	my @columna = @_;
	my $valor = 0;
	my $sumatorio = 0;
	my $elemento; 
	
	foreach $elemento (@columna){
		chomp $elemento;
		$valor = $elemento;
		$sumatorio = $sumatorio + $valor;	
	}
	return $sumatorio;
}