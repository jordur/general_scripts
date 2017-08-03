#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic GATK in target regions
# @Author: JM Rosa
# @Contributors:
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

sub getpos ($$);
sub get_info ($$);
sub get_type ($$);
sub get_genotype ($);
sub analysing_field_8 ($$);

my $F3_var_id;
my $F3_gatk_id;
my $pair_var_id;
my $pair_gatk_id;
my $output_id;


GetOptions (	
				"vf=s"	=> \$F3_var_id,	
				"vp=s"	=> \$pair_var_id,	
				"gf=s"	=> \$F3_gatk_id,	
				"gp=s"	=> \$pair_gatk_id,					
				"o=s"	=> \$output_id,	
);

print "\n\n=================================================================
Collecting Somatic Vars from VarScan and GATK in target regions script v 1.0\n";

if (not defined $F3_var_id or not defined $F3_gatk_id or not defined $pair_var_id or not defined $pair_gatk_id) {
	die "
			Options:
		-vf VarScan F3 file
		-vp VarScan pair file
		-gf GATK F3 file
		-gp GATK pair file
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

my $path = `pwd`;
chomp $path;
 
if (not defined $output_id) {
	$output_id = "varscan_gatk_variants";
}
 
open (LOG, ">$output_id.log") or die "Can´t create log file\n";
print LOG "sg_collecting_varscan_gatk_somatic.pl -vf $F3_var_id -vp $pair_var_id -gf $F3_gatk_id -gp $pair_gatk_id -o $output_id in $path\n";
print LOG "\nStarting process... ".localtime()."\n";

open (OUT, ">$output_id.collected.vcf") or die "Can´t create $output_id.scollected.vcf\n";

open (ERR, ">$output_id.error") or die "Can´t create $output_id.error\n";

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDATA\n";

my %Variants;

open (PGATK, "< $pair_gatk_id") or die "Can´t open $pair_gatk_id\n";

while (my $line = <PGATK>) {
	chomp $line;
	if ($line =~ m/^#Chr/) {
		next;
	}
	else {
		my @data = split (/\s+/, $line);
		my $chr = $data[0];
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;
		if ($chr =~ m/\\0\\0\\/) {
			print "$line\n";
		}
		
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{chr} = $data[0];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{type} = $data[4];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_depth} = $data[5];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_geno} = $data[6];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_freq} = $data[7];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_depth} = $data[8];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_geno} = $data[9];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_freq} = $data[10];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{sb} = $data[11];
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{gatk} = 1;
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{pair} = 1;
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 0;
		$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 1;

	}
} 
close PGATK;

open (F3GATK, "< $F3_gatk_id") or die "Can´t create $F3_gatk_id\n";

while (my $line = <F3GATK>) {
	chomp $line;
	if ($line =~ m/^#Chr/) {
		next;
	}
	else {
		my @data = split (/\s+/, $line);
		my $chr = $data[0];
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;
        if ($chr =~ m/\\0\\0\\/) {
            print "$line\n";
        }
		
		if (defined $Variants{$chr}{$data[1]}{$data[2]}{$data[3]}) {
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 1;
		}
		else {
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{chr} = $data[0];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{type} = $data[4];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_depth} = $data[5];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_geno} = $data[6];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_freq} = $data[7];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_depth} = $data[8];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_geno} = $data[9];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_freq} = $data[10];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{sb} = $data[11];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{gatk} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 0;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{pair} = 0;
		}
	}
} 
close F3GATK;

open (VARSC, "< $pair_var_id") or die "Can´t create $pair_var_id\n";

while (my $line = <VARSC>) {
	chomp $line;
	if ($line =~ m/^#Chr/) {
		next;
	}
	else {
		my @data = split (/\s+/, $line);
		my $chr = $data[0];
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;
        if ($chr =~ m/\\0\\0\\/) {
            print "$line\n";
        }

		if (defined $Variants{$chr}{$data[1]}{$data[2]}{$data[3]}) {
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{pair} = 1;
		}
		else {
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{chr} = $data[0];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{type} = $data[4];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_depth} = $data[5];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_geno} = $data[6];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_freq} = $data[7];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_depth} = $data[8];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_geno} = $data[9];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_freq} = $data[10];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{sb} = $data[11];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{pair} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{gatk} = 0;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 0;
		}
	}
} 
close VARSC;

open (F3VARSC, "< $F3_var_id") or die "Can't open $F3_var_id\n";

while (my $line = <F3VARSC>) {
	chomp $line;
	if ($line =~ m/^#Chr/) {
		next;
	}
	else {
		my @data = split (/\s+/, $line);
		my $chr = $data[0];
		$chr =~ s/chr//;
		$chr =~ s/X/23/;
		$chr =~ s/Y/24/;
		$chr =~ s/M/25/;
        if ($chr =~ m/\\0\\0\\/) {
            print "$line\n";
        }

		if (defined $Variants{$chr}{$data[1]}{$data[2]}{$data[3]}) {
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 1;
		}
		else {	
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{chr} = $data[0];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{type} = $data[4];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_depth} = $data[5];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_geno} = $data[6];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{t_freq} = $data[7];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_depth} = $data[8];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_geno} = $data[9];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{n_freq} = $data[10];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{sb} = $data[11];
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{varscan} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{F3} = 1;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{gatk} = 0;
			$Variants{$chr}{$data[1]}{$data[2]}{$data[3]}{pair} = 0;
		}
	}
} 
close F3VARSC;

foreach my $chr (sort {$a <=> $b} keys %Variants) {
	foreach my $position (sort {$Variants{$chr}{$a} <=>$Variants{$chr}{$b}} keys %{$Variants{$chr}}) {
		foreach my $ref (keys %{$Variants{$chr}{$position}}) {
			foreach my $var (keys %{$Variants{$chr}{$position}{$ref}}) {
				print OUT "$Variants{$chr}{$position}{$ref}{$var}{chr}\t$position\t.\t$ref\t$var\t1000\t$Variants{$chr}{$position}{$ref}{$var}{type}\tTD=$Variants{$chr}{$position}{$ref}{$var}{t_depth};ND=$Variants{$chr}{$position}{$ref}{$var}{n_depth};TG=$Variants{$chr}{$position}{$ref}{$var}{t_geno};NG=$Variants{$chr}{$position}{$ref}{$var}{n_geno};TF=$Variants{$chr}{$position}{$ref}{$var}{t_freq};NF=$Variants{$chr}{$position}{$ref}{$var}{n_freq};SB=$Variants{$chr}{$position}{$ref}{$var}{sb}\tGT:GK:VS:PA:F3\t0/1:$Variants{$chr}{$position}{$ref}{$var}{gatk}:$Variants{$chr}{$position}{$ref}{$var}{varscan}:$Variants{$chr}{$position}{$ref}{$var}{pair}:$Variants{$chr}{$position}{$ref}{$var}{F3}\n";
			}
		}
	}
}

close OUT;
close ERR;

print LOG "\nProcess finished... ".localtime()."\n";
close LOG;

print  "\nProcess finished... ".localtime()."\n";
exit;
