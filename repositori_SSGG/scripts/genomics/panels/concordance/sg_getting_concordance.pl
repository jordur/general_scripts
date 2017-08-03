#!/usr/bin/perl -w
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script for obtaining concordance with HapMap described variants
# @Author: JM Rosa & Arbol
# @Contributors:
########################################################

use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Copy;
use Statistics::Descriptive qw(:all);
use Data::Dumper;
use Tie::File;
use Parallel::ForkManager;
use Scalar::Util;
use analysis_info;

sub getting_common ($$);
sub getting_cov ($$);
sub comparison ($$$$);
sub getting_type($);
sub getting_number($$);
sub getting_novel ($$);

my $coverage_id;
my $vcf_id;
my $xls_id;
my $output_id;
my $threads = 2;
my $sample;
my $type = "both";
my $sep_id = '|';

GetOptions (	
				"c=s"	=> \$coverage_id,
				"d=s"	=> \$sep_id,
				"v=s"	=> \$vcf_id,
				"x=s"	=> \$xls_id,
				"s=s"	=> \$sample,
				"t=s"	=> \$type,
				"o=s"	=> \$output_id,
				"n=i"	=> \$threads,
);

print "\n\n=================================================================
Script for obtaining variants concordance against HapMap controlled cell lines v 1.1\n";

if (!$coverage_id || !$vcf_id || !$xls_id || !$sample) {
	die "
			Options:
		-c Coverage file name  
		-v VCF files name prefix for HapMap variants 
		-x CSV results file name
		-s Sample name that will get compared against HapMap cell line
		-t Type of variants to consider (known, unknown or both, default = both)
		-o Output name path (optional)
		-n Number of threads (optional, default 2)
		-d Separator (optional, default = '|')
\n".localtime()."\n=================================================================\n\n";
}

print "\n Starting process... ".localtime()."\n\n";

if (not defined $output_id) {
	my @name = split (/\./, $coverage_id);
	$output_id = $name[0];
}

my $pm = new Parallel::ForkManager($threads); # Run in parallel with 8 threads max

$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; } );
$pm->run_on_wait( sub {}, 0.5 );

my $config = analysis_info->get_config(type => $type, psv_file => $xls_id);

foreach my $chrom (1..23) {
	my $pid = $pm->start($chrom) and next;
	if ($chrom == 23) { $chrom = 'X';}
	print "Starting analysis for chr$chrom...".localtime()."\n";
	my %Known_snps;
	getting_common (\$chrom, \%Known_snps);
	getting_novel (\$chrom, \%Known_snps);
	$pm->finish;
}
$pm->wait_all_children;

sub getting_novel ($$) {
	my ($chr_id, $known_snps) = @_;
	my @SNVs = `grep \"${sep_id}chr$${chr_id}${sep_id}" $xls_id`;
	open (OUT, "> $output_id\_novel.$$chr_id") or die "Can't create $output_id\_novel.$$chr_id\n";
	foreach my $var (@SNVs) {
		chomp $var;
		my @data = split (/[$sep_id]/, $var);
		
		# Get values from variation in csv file
		my $chr = $data[$config->{chr}];
		my $pos = $data[$config->{position}];
		my $ref = $data[$config->{ref_allele}];
		my $alt = $data[$config->{var_allele}];
		my $zygo = $data[$config->{samples}->{$sample}->{Genotype}];
		my $id = $data[$config->{existing_col}];
		my $depth = $data[$config->{samples}->{$sample}->{Depth}];
		my $MQ = $data[$config->{samples}->{$sample}->{Freq}];
		
		my $value = 0;
		if (defined $$known_snps{$chr}{$pos}{allele1}) {
			#print "$$known_snps{$chr}{$pos}{allele1}\n";
			if ($$known_snps{$$chr_id}{$pos}{allele1} eq $ref or $$known_snps{$$chr_id}{$pos}{allele1} eq $alt) {
				$value++;
			}
			if ($$known_snps{$$chr_id}{$pos}{allele2} eq $ref or $$known_snps{$$chr_id}{$pos}{allele2} eq $alt) {
				$value++;
			}
			
			if ($value > 1) {
				next;
			}
			else {
				print OUT "$var\n";
			}
		}
		else {
			print OUT "$var\n";
		}
	}
	close OUT;
}

sub getting_common ($$) {
	my ($chr_id, $known_snps) = @_;
	open (VCF, "< $vcf_id.chr$$chr_id.snv") or die "Can't open $vcf_id.chr$$chr_id.snv\n";
	open (OUT, "> $output_id\_concordance.$$chr_id") or die "Can't create $output_id.$$chr_id\n";
	while (my $line = <VCF>) {
		my ($type, $concordante, $freq, $id, $ref, $var, $cov, $concordance);
		chomp $line;
		# Do not consider the header
		if ($line =~ "^#") {
			next;
		}
		else {
			my @data = split (/\s+/, $line);

            my $coverage = "NA";
            if (defined $data[7]) {
                $coverage = $data[7];
            }

			my @alleles = split (/\//, $data[6]);
			if ($$chr_id eq 'X' or $$chr_id eq 'M') {
				if (scalar (@alleles) == 1) {
					$alleles[1] = $alleles[0];
				}
			}
			
			if ($alleles[0] eq $data[2] && $alleles[1] eq $data[2]) {
				$type = "Homo_ref";
			}
			elsif ($alleles[0] eq $data[3] && $alleles[1] eq $data[3]) {
				$type = "Homo_var";
			}
			elsif ($alleles[0] eq $data[2] && $alleles[1] eq $data[3]) {
				$type = "Hetero";
			}
			elsif ($alleles[0] eq $data[3] && $alleles[1] eq $data[2]) {
				$type = "Hetero";
			}
			$id = $data[5];
			$freq = $data[4];
			$ref = $data[2];
			$var = $data[3];
			
			$$known_snps{$$chr_id}{$data[1]}{allele1} = $ref;
			$$known_snps{$$chr_id}{$data[1]}{allele2} = $var;
				
			$cov = (&getting_cov($chr_id, $data[1]) || '0');
			if ($cov == 0) {
				$concordance = "No_call\t$type\t";
				print OUT "$$chr_id\t$data[1]\t$ref\t$var\t$id\t$freq\t$data[6]\t$concordance\t$coverage\t$cov\n";
			}
			else {
				my $number = (&getting_number($chr_id, \@data) || '0');
				if ($number == 0) {
					my $i = 0;
					$concordance = &comparison($chr_id, \@data, $type, $i);
					print OUT "$$chr_id\t$data[1]\t$ref\t$var\t$id\t$freq\t$data[6]\t$concordance\t$coverage\t$cov\n";
				}
				else {
					for (my $i = 0; $i < $number; $i++) {
						$concordance = &comparison($chr_id, \@data, $type, $i);
						print OUT "$$chr_id\t$data[1]\t$ref\t$var\t$id\t$freq\t$data[6]\t$concordance\t$coverage\t$cov\n";
					}
				}
			}
		}
	}
	close OUT;
}


sub comparison ($$$$){
	my ($chr_id, $data, $type, $id) = @_;
	my $pos = $$data[1];
	my @col = ('$' . ($config->{chr} + 1),'$' . ($config->{position} + 1));
	my @var = `awk -F\'$sep_id\' \'$col[0] == \"chr$${chr_id}\" && $col[1] == $pos\' $xls_id`;
	my $value = 0;
	my $my_type;

	# If no variant has been found for the given position, it is considered Homo_ref (since the coverage for the position was already checked!!)
	if (!@var) {
		$my_type = "Homo_ref";
		if ($my_type eq $type) {
			return "concordant\tHomo_ref\t-";
		}
		else {
			if ($type eq "Homo_var") {
				return "discordant\tHomo_var\tHomo_ref";
			}
			else {
				return "discordant\tHetero\tHomo_ref";
			}
		}
	}
	
	else {
		$my_type = &getting_type($var[$id]);		
		if ($my_type eq $type) {
			if ($type eq "Homo_var") {
				return "concordant\tHomo_var\t-";
			}
			elsif ($type eq "Homo_ref") {
				return "concordant\tHomo_ref\t-";
			}
			else {
				return "concordant\tHetero\t-";
			}
		}
		else {
			if ($type eq "Homo_var" && $my_type eq "Hetero") {
				return "discordant\tHomo_var\tHetero";
			}
			elsif ($type eq "Homo_var" && $my_type eq "Homo_ref") {
				return "discordant\tHetero\tHomo_var";
			}
			elsif ($type eq "Hetero" && $my_type eq "Homo_var") {
				return "discordant\tHetero\tHomo_var";
			}
			elsif ($type eq "Hetero" && $my_type eq "Homo_ref") {
				return "discordant\tHetero\tHomo_ref";
			}
			elsif ($type eq "Homo_ref" && $my_type eq "Hetero") {
				return "discordant\tHomo_ref\tHetero";
			}
			elsif ($type eq "Homo_ref" && $my_type eq "Homo_var") {
				return "discordant\tHomo_ref\tHomo_var";
			}
			else {
				return "Unknown_discordance\t$type\t$my_type";
			}
		}
	}
}

sub getting_number($$) {
	my ($chr_id, $data) = @_;
	my $pos = $$data[1];
	my @col = ('$' . ($config->{chr} + 1),'$' . ($config->{position} + 1));
	my @var = `awk -F\'$sep_id\' \'$col[0] == \"chr$${chr_id}\" && $col[1] == $pos\' $xls_id`;
	my $number = scalar @var;
	return $number;
}

sub getting_type($) {
	my $var = shift;
	my @data = split (/[$sep_id]/, $var);
	my $type;
	if (Scalar::Util::looks_like_number($data[$config->{samples}->{$sample}->{Freq}])){
		if ($data[$config->{samples}->{$sample}->{Freq}] < 0.15) {
			$type = "Homo_ref";
		}
		elsif ($data[$config->{samples}->{$sample}->{Freq}] > 0.36 and $data[$config->{samples}->{$sample}->{Freq}] < 0.65) {
			$type = "Hetero";
		}
		elsif ($data[$config->{samples}->{$sample}->{Freq}] > 0.85) {
			$type = "Homo_var";
		}
		else {
			$type = "Uncertain";
		}
	} else {
		$type = "Uncertain";
	}
	return $type;
}

sub getting_cov ($$) {
	my ($chr, $pos) = @_;
	my $col = '$4';
	my $cov = `grep -m 1 "chr$$chr\[\[:space:\]\]$pos\[\[:space:\]\]" $coverage_id.chr$$chr.pileup | awk \'{print $col}\'`;
	chomp $cov;
	#print "$$chr\t$pos\t$cov\n";
	return $cov;
}

print "Process finished... ".localtime()."\n\n";

exit;
