#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (sum);
use File::Basename;

### Input Data ###
#CSV with gene list
my $input_csv_gene="/home/gmarco/Onco2/oncogeneprofile.csv";

#Gene BED file directory

my $input_bed_not_covered = "";
my $input_gene_directory= "/home/gmarco/Onco2/ONCOGENEPROFILE_min0/genes";
my $output_bed_gene_directory= "/home/gmarco/Onco2/ONCOGENEPROFILE_min0/genes/bed";
my $input_bed_gene_directory = $output_bed_gene_directory;

## Get gene list
my @gene_list = &get_gene_list($input_csv_gene);
#&convert_bed($input_gene_directory, $output_bed_gene_directory, @gene_list);

main();

sub main {
	my $gene;
	my ($chr,$input_totals_coords,$input_exon_coords,$input_splicing_coords,$input_strong_splicing_coords,$input_UTR_coords,$input_vars_coords);
	
	## Open output_file
	open OUTPUT, ">/home/gmarco/Onco2/output_guille";
	print OUTPUT "HGNC\tCHR\tNC_EXON\tNC_SPLICING\tNC_STRONG_SPLICING\tNC_UTR\tNC_VARS\n";
	
	foreach $gene (@gene_list){
		chomp $gene;
 		
		$input_totals_coords = $input_bed_gene_directory."/$gene"."_totals_coords.bed";
		$input_exon_coords = $input_bed_gene_directory."/$gene"."_exon_coords.bed";
		$input_splicing_coords = $input_bed_gene_directory."/$gene"."_splicing_coords.bed";
		$input_strong_splicing_coords = $input_bed_gene_directory."/$gene"."_strong_splicing_coords.bed";
		$input_UTR_coords = $input_bed_gene_directory."/$gene"."_UTR_coords.bed";
		$input_vars_coords = $input_bed_gene_directory."/$gene"."_vars_coords.bed";
		
		$chr = (split ("\n",`awk '{print \$1}' $input_totals_coords`))[0];
		#my @columna_exon = `subtractBed -a $input_totals_coords -b $input_exon_coords | awk '{print \$3-\$2+1}'`;
		
		my $substract_command = "subtractBed -a $input_totals_coords -b ";
		my $awk_command = " | awk '{print \$3-\$2+1}'";
		my $awk = "awk '{print \$3-\$2+1}'";
		
		my $total_gen = &col_sum(`$awk $input_totals_coords`);
		my $total_exon_coords = &col_sum(`$substract_command $input_exon_coords $awk_command`);
		my $total_splicing_coords = &col_sum(`$substract_command $input_splicing_coords $awk_command`);
		my $total_strong_splicing_coords = &col_sum(`$substract_command $input_strong_splicing_coords $awk_command`); 
		my $total_UTR_coords = &col_sum(`$substract_command $input_UTR_coords $awk_command`);
		my $perc_vars;
		
		if (-z "$input_vars_coords") {
			my $total_vars_coords = 0;
			$perc_vars = 0;
		}
		else {
			my $total_vars_coords =  &col_sum(`$substract_command $input_vars_coords $awk_command`);
			$perc_vars = &percentage_coverage($total_gen, $total_vars_coords);
		}
		
		my $perc_exon = &percentage_coverage($total_gen, $total_exon_coords);
		my $perc_splicing = &percentage_coverage($total_gen, $total_splicing_coords);
		my $perc_strong_splicing = &percentage_coverage($total_gen, $total_strong_splicing_coords);
		my $perc_UTR = &percentage_coverage($total_gen, $total_UTR_coords);
		
		print OUTPUT "$gene\t$chr\t$perc_exon\t$perc_strong_splicing\t$perc_strong_splicing\t$perc_UTR\t$perc_vars\n";

	}
	close OUTPUT; 
}

sub percentage_coverage {
	my ($total_gen, $total_columna) = @_;
	my $porc = 100*(1 - ($total_columna/$total_gen));
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

sub convert_bed{
	open CONVERT_COMMANDS, ">/home/gmarco/Onco2/convert_commands.sh";
	my ($input_gene_directory,$output_bed_gene_directory, @gene_list) = @_;

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