#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic GATK in target regions
# @Author: JM Rosa
# @Contributors: Sheila
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
use genomics::genotype;

sub getpos ($$);
sub get_info ($$$);
sub analysing_field_8 ($$);

my $vcf_id;
my $output_id;
my $sample_id;

GetOptions (	
				"v=s"	=> \$vcf_id,	
				"o=s"	=> \$output_id,	
				"s=s"	=> \$sample_id,	
);

print "\n\n=================================================================
Getting Vars from GATK in target regions script v 1.0\n";

if (not defined $vcf_id or not defined $sample_id) {
	die "
			Options:
		-v VCF file
		-s Sample id (comma separated)
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;
 
if (not defined $output_id) {
	my $name = (split(/\./, $vcf_id))[0];
	$output_id = $name."_gatk_parsed";
}
 
open (LOG, ">$output_id.log") or die "Can't create log file\n";
print LOG "sg_parsing_gatk.pl -v $vcf_id -s $sample_id -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

open (VCF, "< $vcf_id") or die "Can't open $vcf_id\n";
my %columns;
open (OUT, ">$output_id.snvs") or die "Can't create $output_id.snvs\n";
open (ERR, ">$output_id.error") or die "Can't create $output_id.error\n";

my @samples = split (/,/, $sample_id);
print OUT "#Chr\tPosition\tRef_Allele\tVar_Allele\tType\t";
foreach my $sample (@samples) {
	print OUT "$sample\_depth\t$sample\_geno\t$sample\_freq\t";
}
print OUT "SB\n";

while (my $line = <VCF>) {
	chomp $line;
	if ($line =~ m/^##/) { next;}
	elsif ($line =~ m/^#CHROM/) {
		foreach my $sample (@samples) {
			$columns{$sample} = &getpos ($sample, $line);
			if ($columns{$sample} eq ""){
				die("ERROR: Sample $sample not found in vcf!!\n");
			}
		}
	}
	else {
		my @data = split (/\s+/, $line);
		my $type = '-';
		my @vars = split (/\,/, $data[4]);
		my $n = 1;
		my $colinfo= $data[8];
		for (my $var=0; $var<=$#vars ; $var++)
		{
			print OUT "$data[0]\t$data[1]\t$data[3]\t$vars[$var]\t$type\t";
			my ($depth, $info, $freq, $geno);
			my $allele=$var+1;
			foreach my $sample (@samples) {
				if (($data[$columns{$sample}] ne './.') && ($data[$columns{$sample}] ne '.') && ($data[$columns{$sample}] ne '0/0')){
					$depth = &get_info ($data[$columns{$sample}], 'depth', $n, $colinfo);
					$info = &get_info ($data[$columns{$sample}], 'freq', $n, $colinfo);
					$freq = (split(/\,/, $info))[1];
					$geno = (split(/\,/, $info))[0];	
				}
				else {
					$depth = '-';
					$info = '-,-';
					$freq = (split(/\,/, $info))[1];
					$geno = (split(/\,/, $info))[0];	
				}
				print OUT "$depth\t$geno\t$freq\t";
			}
			# Since no information on the SB parameter is available, the FS parameter
			# will be instead considered for the strand bias estimation. FS is the phred
			# scale p-value of the Fisher's exact test of the contingency table between 
			# reference alleles count in the + and - strand and variant alleles count
			# in the + and - strand:
			my $sb = &analysing_field_8 ($data[7], 'FS=');
			print OUT "$sb\n";
			$n++;
		}	
	}
} 

close VCF;
close OUT;
close ERR;

sub analysing_field_8 ($$){
	my ($line, $match) = @_;
	my @fields =  split (/;/, $line);
	my $value = '-';
	foreach my $field (@fields) {
		if ($field =~ m/$match/) {
			$field =~ s/$match//;
			$value = $field;
		}
	}
	return $value;
}

sub get_info ($$$) {
	my ($line, $match, $n, $format) = @_;
	my @fields =  split (/:/, $line);
	my $var;
	my $ref;
	my @split_format= split (/:/,$format);
	for (my $i=0;$i<=$#split_format;$i++)
	{
		if ($split_format[$i] eq "AD"){
			$ref = (split (/\,/, $fields[$i]))[0];
			$var = (split (/\,/, $fields[$i]))[1];
		}	
	}
	my $value;

	if ($match eq 'depth') {
		$value = $var + $ref;
	}
	elsif ($match eq 'freq') {
		my $depth = $var + $ref;
		if ($depth == 0) {
			$value = '-,-';
		}
		else {
			my $freq = $var / $depth;
			$value = genotype->get_genotype ($depth,$freq,$ref,$var); 
		}
	}
	return $value;
}

sub getpos ($$) {
	my ($id, $line) = @_;
	my @data = split(/\s+/,$line);
	for (my $i = 0; $i <= $#data; $i++){
		my $sample = substr($data[$i],0,length($id));
		if ($sample eq $id) {
			return $i;
			last;
		}
	}
}

print LOG "\nProcess finished... ".localtime()."\n";
close LOG;

print  "\nProcess finished... ".localtime()."\n";
exit;
