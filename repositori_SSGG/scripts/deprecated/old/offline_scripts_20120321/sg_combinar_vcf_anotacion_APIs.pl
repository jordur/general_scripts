#!/usr/bin/perl

#==============================#
# Combining anot and vcf files #
#==============================#


use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Copy;
#use Parallel::ForkManager; 
#use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;

sub es_indel ($);
sub compare_alleles ($$$$$);
sub homoheterozigoto ($);
sub longitud_variacion ($);

my ($vcf_id, $anot_id) = @ARGV;

open (VCF, "< $vcf_id") or die "Can´t open $vcf_id\n";

my @VCF = <VCF>;

close VCF;

open (ANOT, "< $anot_id") or die "Can´t open $anot_id\n";

while (my $line = <ANOT>) {
	chomp $line;
	my @anot = split (/\t/, $line);
	$anot[0] =~ s/chr//;
	my $type = es_indel($anot[2]);

	my $position;
	
	if ($type eq 'indel') {
		$position = $anot[1] - 1;
	}
	else {
		$position = $anot[1];
	}
	foreach my $variant (sort @VCF) {
		chomp $variant;
		my @vcf = split (/\t/, $variant);
		
		$vcf[0] =~ s/chr//;
		$vcf[7] =~ s/;DS;/;/;
		my @col8 = split (";",$vcf[7]);
		$col8[1] =~ s/AF=//;
		$col8[3] =~ s/DP=//;
		$col8[7] =~ s/MQ=//;
		
		my @col10 = split (":",$vcf[9]);
		
		if (($vcf[0] eq $anot[0]) && ($vcf[1] == $position)) {
			my $index = compare_alleles($anot[2], $vcf[3], $vcf[4], $line, $variant);
		
			if ($index == 1){

				$anot[13]=~s/\./,/;
				$col10[3]=~s/\./,/;
				$col8[7]=~s/\./,/;	
				
				my $status = homoheterozigoto($col8[1]);
				my $length = longitud_variacion ($anot[2]);
				
				print $anot[6],"\t",$vcf[0],"\t",$vcf[1],"\t",$vcf[3],"\t", $vcf[4],"\t",$type,"\t",$length,"\t",$status,"\t",$anot[12],"\t",$anot[4],"\t",$anot[5],"\t",$anot[7],"\t",$anot[8],"\t",$anot[9],"\t",$anot[10],"\t",$anot[11],"\t",$anot[13],"\t",$anot[14],"\t",$anot[15],"\t",$col8[3],"\t", $col10[2],"\t",$col8[7],"\t",$anot[16],"\t",$anot[17],"\n";

				last;
			}
			else {
				next;
			}
		}
		else {
			next;
		}
	}
}

close ANOT;

sub compare_alleles ($$$$$) {
	my ($genotype, $ref_vcf, $var_vcf, $line, $variant) = @_;
	my $ref_anot = (split(/\//, $genotype))[0];
	my $var_anot = (split(/\//, $genotype))[1];
	
	if ($ref_anot eq $ref_vcf && $var_vcf eq $var_anot) {
		return 1;
	}
	
	else {
		my $trans_ref = $ref_vcf;
		$trans_ref =~ s/$ref_anot//;
		my $trans_var = $trans_ref.$var_anot;
		$trans_var =~ s/-//;
		
		if ($var_vcf eq $trans_var) {
			return 1;
		}
		else {
			return 0;
		}
	}

}

sub longitud_variacion ($) {
	my ($string) = @_;
	my @split_variation = split ("/",$string);
	my $lalong = abs(length($split_variation[0])-length($split_variation[1]));

	if ($split_variation[0] eq "-") {
		return (length ($split_variation[1]));
	}	
	elsif ($split_variation[1] eq "-") {
		return (length ($split_variation[0]));
	}

	elsif (length ($split_variation[1]) eq length ($split_variation[0])) {
		return "-";
	}
	elsif (($split_variation[0] ne "-") && ($split_variation[1] ne "-") && ($lalong > 0) ) {
		return $lalong;
    }
}

sub homoheterozigoto ($) {
	my ($status) = @_;
	if ($status eq "0.50") {
		return "heterozygous";
	}
	else {
		return "homozygous";
	}
}

sub es_indel ($) {
	my ($genotype) = @_;
    my @split_variation = split ("/",$genotype);
    if (($split_variation[0] eq "-") || ($split_variation[1] eq "-") || (length($split_variation[0]) != length ($split_variation[1]))) {		
		return "indel";
	}
	else
	{
		return "SNV";
	}
	
}

exit;