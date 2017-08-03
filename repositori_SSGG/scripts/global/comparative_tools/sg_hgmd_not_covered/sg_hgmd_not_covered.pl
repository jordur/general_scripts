#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: HGMD(Biobase) non covered design zones
# @Author: Guillermo Marco Puche
# @Contributors: 
########################################################

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Getopt::Std;

main();

sub main {
	
	my $masterbiobed = "./masterbioBed_sorted.bed";
	my $bed_base_input;
	my %options=();
	
	getopts("r:i:o:s:", \%options);
	
	print "\n=================================================================
	HGMD NON COVERED REGIONS\n";
	#print "\t",`date`;
	if (not defined $options{r} or not defined $options{i} or not defined $options{o}or  not defined $options{s}){
		print "\t!!! YOU MUST SPECIFY ALL INPUT FLAGS !!!\n\n";
		print "-r *.bed target_regions file (BED) USING ABSOLUTE PATH !!! ie: /home/data/target_regions.bed\n";
		print "-i BAM alignment file USING ABSOLUTE PATH !!! ie: /home/mapping/file.bam\n";
		print "-o Output path USING ABSOLUTE PATH !!! ie /home/resulsts/output/\n";
		print "-s Sample name !!\n";
		print "=================================================================\n\n";
		exit;
	}
	print "=================================================================\n";

	#Check input region file exists
	if (not -e $options{r}){
		die "ERROR: Couldn't read BED target regions file.\n";
	}
	
	#Check input BAM align file exists
	if (not -e $options{i}){
		die "ERROR: Couldn't read BAM align file.\n";
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

	#$bed_base_input = basename($options{r},".bed");
	$bed_base_input = basename($options{s});
	
	#Convert BAM input alignment to BED.
	#bamToBed va fuera ya que usamos el bed en el primer script con samtools depth
	print "Converting BAM alignment to BED and calculating non covered design regions..\n";
	system("bamToBed -ed -i $options{i} | subtractBed -a $options{r} -b stdin > $output_dir/$options{s}_not_covered.bed");
	print "DONE ! Not covered design regions regions obtained: $output_dir/$options{s}_not_covered.bed\n";
	print "=================================================================\n";
	
	print "Obtaining (HGMD) information for non covered regions...\n";

	#Obtain HGMD information for non covered regions...
	print `intersectBed -a $output_dir/$options{s}_not_covered.bed -b $masterbiobed -wa -wb | awk 'BEGIN { OFS = "\\t";} {print \$1,\$2,\$3,\$7}' > $output_dir/$options{s}_not_covered_hgmd_mutations.bed`;
	
	print "DONE ! HGMD information for non covered design regions regions: $output_dir/$options{s}_not_covered_hgmd_mutations.bed\n";
	print "=================================================================\n";
	print "Obtaining number of bases for each non covered region...\n";
	
	#Obtain number of bases in column format for non_covered regions.
	system("awk '{print \$3-\$2+1}' $output_dir/$options{s}_not_covered.bed > $output_dir/$options{s}_not_covered_bases_column.bed");
	
	print "DONE ! Number of bases for each non covered region: $output_dir/$options{s}_not_covered_bases_column.bed\n";
	print "=================================================================\n";
	
	#Merge non_covered regions with bases count for each line
	print `paste -d'\\t' $output_dir/$options{s}_not_covered.bed $output_dir/$options{s}_not_covered_bases_column.bed > $output_dir/$options{s}_not_covered_merged_bases.bed`;
	
	my $base_merged_bed = "$output_dir/$options{s}_not_covered_merged_bases.bed";
	my $hgmd_bed = "$output_dir/$options{s}_not_covered_hgmd_mutations.bed";
	my $final = "$output_dir/$options{s}"."_final.bed";
	
	print "Merging final resulsts file...\n";
	final_merge($base_merged_bed, $hgmd_bed, $final);
	print "DONE ! Final resulsts file: $output_dir/$options{s}"."_final.bed\n";
	print "=================================================================\n";
	print "HGMD NON COVERED REGIONS SCRIPT SUCCEED !\n";
}

sub final_merge {
	my ($base_merged_bed, $hgmd_bed, $final) = @_;
	
	open BASE, "$base_merged_bed";
	open FINAL, ">$final";
	
	while (<BASE>){
		chomp $_;
		my $found = 0;
		my @split_base = split("\t", $_);

		open HGMD, "$hgmd_bed";
		
		while (<HGMD>){
			chomp $_;
			my @split_hgmd = split("\t", $_);
			if ($split_base[0] eq $split_hgmd[0] && $split_base[1] == $split_hgmd[1] && $split_base[2] == $split_hgmd[2]){
				$found = 1;
				print FINAL "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\t$split_hgmd[3]\n";
			} 	
		}
		close HGMD;
		
		if ($found == 0) {
			print FINAL "$split_base[0]\t$split_base[1]\t$split_base[2]\t$split_base[3]\n";
		}	
	}
	close HGMD;
	close BASE;
	close FINAL;
}

#paste -d"\t" /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/test/panel.bed /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/test/bases_panel.bed

#sub get_bed_noncovered(){
#	my $bed_non_covered = $_[0];
#	my ($chr, $ii, $fi, $start, $end);
#	
#	#We open Non covered BED intervals
#	open BED, "$bed_non_covered";
#	while (<BED>){
#		chomp $_; 
#		my @line = split ("\t", $_);
#		$chr = $line[0]; $ii = $line[1]; $fi = $line[2];
#		#&get_biobase_info($chr, $ii, $fi);
#		system("intersectBed -a $bed_non_covered -b biobed_test.bed -wa -wb | awk '{print $1,$2,$3,$7}'");
#		
#	}
#}

#sub get_bio_file_handle {
#	my $chr = shift;
#	my $bio_file_handle = new FileHandle;
#	my $file = "/home/likewise-open/SGNET/gmarco/biobase/bioBaseSortUniques_chr$chr.txt";
#	$bio_file_handle->open($file) or die("ERROR: Could not read from input file: $file\n");
#	return $bio_file_handle;
#}


#sub convert_biobase_bed () {
#	
#    my ($chr) = $_[0];
#    my $bio_file_handle = &get_bio_file_handle($chr);
#    my $desc = "-\t-";
#    open (OUTPUT, ">/home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/bioBed_chr$chr.bed");
#    while (my $line = <$bio_file_handle>) {
#        chomp $line;
#        
#        if ($line =~ /^Type/) {next;}
#        
#        my @data = split (/\t/, $line);
#        
#        my ($chr, $start, $end) = split (/[:-]/, $data[4]);
#        
#        if (!$end or $end eq '' or $end eq '+') { $end = $start + 1;}
#        
#        $desc = $chr."\t".$start."\t".$end."\t".$data[1]."|".$data[2]."|".$data[6]."|".$data[9]."|".$data[7]."|";
#        #$desc =~ s/\s+/_/g;
#		print OUTPUT $desc,"\n";
#    }
#    system("sortBed -i /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/bioBed_chr$chr.bed > /home/likewise-open/SGNET/gmarco/sg_hgnc_no_cover_testing/bioBed/sorted/bioBed_sorted_chr$chr.bed");
#    close OUTPUT;
#    
#   	#Sort bioBed file
#   	
#    return 1;
#    #return $desc;
#}
