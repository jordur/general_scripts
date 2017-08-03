#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Script to Parsing Somatic GATK in target regions
# @Author: JM Rosa
# @Contributors:
########################################################

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use genomics::genotype;

open (GFF,"<",$ARGV[0]);

print "#Chr\tPosition\tRef_Allele\tVar_Allele\t$ARGV[1]_depth\t$ARGV[1]_geno\t$ARGV[1]_freq\tReads+\tReads-\n";
while (my $lines=<GFF>)
{
	chomp ($lines);
	my @spltlinea=split ("\t",$lines);
	my $char = substr ($spltlinea[16],0,1);
	my @splitAlleles=split ("/",$spltlinea[9]);
	if ($spltlinea[2] eq "deletion")
	{
		# Positions from BioScope gff have an offset of 1
		print "$spltlinea[0]\t",$spltlinea[3]-1,"\t$spltlinea[16]\t$char\t$spltlinea[12]\t",&zygosity($spltlinea[12],$spltlinea[14],$spltlinea[16],$char),"\t",&calculate_freq($spltlinea[14]),"\t",&strandBias($spltlinea[15]),"\n";
	}
	elsif (($spltlinea[2] eq "insertion_site"))
	{
		print "$spltlinea[0]\t$spltlinea[3]\t$char\t",$char.$splitAlleles[1],"\t$spltlinea[12]\t",&zygosity($spltlinea[12],$spltlinea[14],$char,$char.$splitAlleles[1]),"\t",&calculate_freq($spltlinea[14]),"\t",&strandBias($spltlinea[15]),"\n";
	}

}

sub calculate_freq
{
	my $zygosity_score = shift;
	
	# BioScope has a zygosity-score with value 0 in case of homozygotes for the variant, and value 1 in case of heterozygotes.
	# For normalizing with the SG frequency scale, following transformation is needed: SG_freq = -0.5 * zygosity_score + 1
	
	return -0.5 * $zygosity_score + 1;
}

sub zygosity
{
	my ($depth,$zygosity_score,$ref,$var) = @_;
	
	my $freq = -0.5 * $zygosity_score + 1;
	
	# Different criteria have been stablished for indels and SNVs (summarized in genotype.pm library)
	my @tmp = split(",", genotype->get_genotype ($depth,$freq,$ref,$var));
	my $value = $tmp[0];
	return $value;
}

sub hemihomo
{
	my $valor = shift;
	if ($valor eq "HOMOZYGOUS")
	{
		return "P_Homo";
	}
	else
	{
		return "P_Hetero";
	}
}

sub strandBias
{
	my $cadenas= shift;
	my $negativos=0;
	my $positivos=0;
	my @splt_cadenas=();
	@splt_cadenas= split (",",$cadenas);
	for (my $a=0;$a<=$#splt_cadenas;$a++)
	{
			if ($splt_cadenas[$a] eq "-")
			{
				$negativos=$negativos+1;
			}
			elsif ($splt_cadenas[$a] eq "+")
			{
				$positivos=$positivos+1;
			}
	}
	return $positivos."\t".$negativos;
}
