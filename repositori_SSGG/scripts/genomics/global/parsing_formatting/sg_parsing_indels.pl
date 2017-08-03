#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script for getting and collecting BioScope indels from different samples in target regions
# @Author: Arbol
# @Contributors: JM Rosa
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
use genomics::variants;

my @vcf_ids;
my $output_id;
my $sample_id;

GetOptions (
				"v=s"	=> \@vcf_ids,
				"o=s"	=> \$output_id,
);

print "\n\n=================================================================
Getting and collecting BioScope indels from different samples in target regions v1.0\n";

if (not defined $vcf_ids[0]) {
	die "
			Options:
		-v VCF files (comma separated list of vcf files, simply append them: -v file1,file2,file3...)
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;

if (not defined $output_id) {
	$output_id = "all_indels_parsed";
}

# Open log, output and error files
open (OUT, ">$output_id.snvs") or die "Can't create $output_id.snvs\n";
open (ERR, ">$output_id.error") or die "Can't create $output_id.error\n";
open (LOG, ">$output_id.log") or die "Can't create log file\n";

print LOG "sg_parsing_gatk.pl -v @vcf_ids -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

# process all vcf files
@vcf_ids = split(/,/,join(',',@vcf_ids));

# First, get all sample names and set file pointers to each file
my @samples; # Array containing the sample names of the files
my @files; # Array containing the handles to the vcf files
my @lines; # Array containing one line of each vcf file
my @variant; # Array containing the corresponding variant
my $SB; #Strand bias value for each indel
my $samples_for_indel; #Number of samples with a given indel

for (my $i=0;$i<=$#vcf_ids;$i++){
	my $vcf_id = $vcf_ids[$i];
	my $header = `head -1 $vcf_id`;
	my @cols = split('\t',$header);
	if (not defined $cols[5]){
		print "WARNING: Vcf file $vcf_id apparently hasn't right format and won't be considered\n";
		push (@samples,"");
	} else {
		my $tmp = substr($cols[4],0,length($cols[4])-6);
		open ($files[$i], "< $vcf_id") or die "ERROR: Can't open $vcf_id\n";
		push (@samples,$tmp);
		my $temp_line;
		do {
			$lines[$i] = readline($files[$i]);
			$temp_line = $lines[$i];
		} while ($temp_line =~ m/^#Chr/);
		$variant[$i] = variants->create_from_tsv($lines[$i]);
	}
}

# Print output header
print OUT "#Chr\tPosition\tRef_Allele\tVar_Allele\tType";
foreach my $sample (@samples) {
	print OUT "\t$sample\_depth\t$sample\_geno\t$sample\_freq";
}
print OUT "\tSB\n";

# Now parse all sample files and output indels
my $first_pos;
my $first_variant;
do {
	$first_variant = variants->get_first_variant(\@variant);
	if (defined($first_variant)){
		my $line = "chr$first_variant->{chr}\t$first_variant->{position}\t$first_variant->{ref_allele}\t$first_variant->{var_allele}\t-";
		$SB = 0;
		$samples_for_indel = 0;
		for (my $i=0;$i<=$#lines;$i++){
			if (defined($variant[$i])){
				if (variants->compare($variant[$i],$first_variant) == 0){
					$line = $line . "\t$variant[$i]->{depth}\t$variant[$i]->{genotype}\t$variant[$i]->{freq}";
					$samples_for_indel ++;
					$SB = $SB + $variant[$i]->{SB};
		
					# Read new line from file and obtain new variant
					$lines[$i] = readline($files[$i]);
					if (defined($lines[$i])){
						$variant[$i] = variants->create_from_tsv($lines[$i]);
					} else {
						$variant[$i] = undef;
					}
				} else {
					$line = $line . "\t-\t-\t-";
				}
			} else {
				$line = $line . "\t-\t-\t-";
			}
		}
		$line = $line ."\t" . ($SB/$samples_for_indel) . "\n";
		print OUT $line;
	}
} while (defined($first_variant));

# Close files
for (my $i=0;$i<=$#vcf_ids;$i++){
	close $files[$i];
}
close OUT;
close ERR;

# Finish log
print LOG "\nProcess finished... ".localtime()."\n";
close LOG;
print  "\nProcess finished... ".localtime()."\n";
exit;