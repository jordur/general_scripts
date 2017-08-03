#!/usr/bin/perl -w
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to compare sample variants against HapMap samples that are in the same VCF file
# @Author: Sheila
# @Contributors: Arbol
# @Date: May 2013
# @Version: 0.2
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Cwd;
use file_handle;

# ---------------------------------
# -- global & default parameters --
# ---------------------------------

my $file_name;
my @split_lines;
my @comparisons;
my @split_pairs;
my @splitSampleValues;
my @splitHapMapValues;
my @header;
my $output;
my $input;
my @sampleNames;
my $hapmapname;
# ---------------------------
# -- functions definitions --
# ---------------------------

sub compare_samples {
    my ($value)= shift;
    
    my (@ctp,@dtp,@tn,@fp,@fn);
    
	open (VCF,"<","$value") or die "ERROR: impossible to open $value!!!\n";
	
	# Obtain data from header
	my $header = 1;
	while (defined(my $lines=<VCF>) and $header){
		chomp($lines);
		@split_lines= split("\t",$lines);
		if ($lines =~ /##/){
			push (@header,$lines);
		}
		if ($lines =~ /#CHROM/){
			$hapmapname=$split_lines[-1];
			for (my $i=9; $i<=$#split_lines-1; $i++){
				push (@comparisons,$i.",".$#split_lines);
				push (@sampleNames,$split_lines[$i]);
				# open output files
				my $file_handle;
				$file_handle = &get_out_file_handle("$output"."CTP_"."$split_lines[$i]".".vcf");
				push(@ctp,$file_handle);
				$file_handle = &get_out_file_handle("$output"."DTP_"."$split_lines[$i]".".vcf");
				push(@dtp,$file_handle);
				$file_handle = &get_out_file_handle("$output"."TN_"."$split_lines[$i]".".vcf");
				push(@tn,$file_handle);
				$file_handle = &get_out_file_handle("$output"."FP_"."$split_lines[$i]".".vcf");
				push(@fp,$file_handle);
				$file_handle = &get_out_file_handle("$output"."FN_"."$split_lines[$i]".".vcf");
				push(@fn,$file_handle);
			}
			$header =0;
		}
	}
	close (VCF);
	open (VCF,"<","$value");
	
	while (my $lines=<VCF>){
	 	chomp($lines);
		@split_lines= split("\t",$lines);
		if (($lines !~ /#CHROM/) && ($lines !~ /##/)){
			for (my $a=0;$a<=$#comparisons;$a++){
				@split_pairs= split(",",$comparisons[$a]);
				@splitSampleValues= split (":",$split_lines[$split_pairs[0]]);
				@splitHapMapValues= split (":",$split_lines[$split_pairs[1]]);
				if ((($splitSampleValues[0] eq $splitHapMapValues[0]) && ($splitHapMapValues[0] eq "0/0"))
					|| (($splitSampleValues[0] eq ".") && ($splitHapMapValues[0] eq "0/0"))  )
				{
					my $fh = $tn[$a];
					print $fh "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
				elsif (($splitSampleValues[0] eq $splitHapMapValues[0]) && ($splitHapMapValues[0] ne "0/0")
					&& ($splitSampleValues[0] ne "0/0") && ($splitHapMapValues[0] ne ".") 
					&& ($splitSampleValues[0] ne "."))
				{
					my $fh = $ctp[$a];
					print $fh "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
				elsif (($splitSampleValues[0] ne $splitHapMapValues[0]) && ($splitHapMapValues[0] ne "0/0")
					&& ($splitSampleValues[0] ne "0/0") && ($splitHapMapValues[0] ne ".") 
					&& ($splitSampleValues[0] ne "."))
				{
					my $fh = $dtp[$a];
					print $fh "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
				elsif ((($splitHapMapValues[0] eq ".") || ($splitHapMapValues[0] eq "0/0")) 
					&& ($splitSampleValues[0] ne ".") &&  ($splitSampleValues[0] ne "0/0"))
				{
					my $fh = $fp[$a];
					print $fh "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
				elsif (($splitHapMapValues[0] ne "0/0") && ($splitHapMapValues[0] ne ".") 
					&& (($splitSampleValues[0] eq "." ) || ($splitSampleValues[0] eq "0/0" )))
				{
					my $fh = $fn[$a];
					print $fh "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
				else {
					print "$split_lines[0]\t$split_lines[1]\t$split_lines[2]\t$split_lines[3]\t$split_lines[4]\t$split_lines[5]\t$split_lines[6]\t$split_lines[7]\t$split_lines[8]\t$split_lines[$split_pairs[0]]\t$split_lines[-1]\n";
				}
			}
		}
	}
	foreach (@ctp) {close($_);}
	foreach (@dtp) {close($_);}
	foreach (@tn) {close($_);}
	foreach (@fp) {close($_);}
	foreach (@fn) {close($_);}
}


# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"i=s"	=> \$input,
				"o=s"	=> \$output,
);

print "\n=================================================================
$basename: Script to compare sample variants from different samples against reference genotype data (usually HapMap lines from 1000Genomes project). Samples must be in the same VCF file, reference data in the last column\n";

if (not defined $input) {
 	die "
			Options:
			-i full path to VCF file that includes all samples and the HapMap cell line data to be compared\n
			-o output path\n
			
			In output path, following files will be created:
			CTP (true positives and concordant in genotype against reference), DTP (true positives but discordant in genotype against reference), FN (false negatives), FP (false positives), TN (true positives)
			
			***NOTE: VCF file needs to have a proper header, the script will look for samples names in the header, at least this line should be present\n
			Ex:	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12144_Replicate1      NA12144_Replicate2      NA12144
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "\nProcess finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	print "Running $basename with parameters i=$input o=$output\n";
	&compare_samples($input)
}

exit;
