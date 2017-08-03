#!/usr/bin/perl -w
#####################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to transform pseudoVCF files into real VCF files
# @Author: Sheila
# @Date: May 2013
# @Version: v0.1
#####################################################

use strict;

### Variable definition
my $lines;
my @split_lines;
my @split_col8;
my $isnotheader="";
my $genotype;
my $allelefreq;
my $depth;
my @sampleInfo;
my $ecRef;
my $ecAlt;
my $samplelist;
my @split_samplelist;
my $valCol;
my $value;
open (PSEUDOVCF,"<",$ARGV[0]);
open (VCF,">>",$ARGV[1]);
$samplelist=$ARGV[2];
while ($lines=<PSEUDOVCF>)
{
	chomp ($lines);
	@split_lines=split("\t",$lines);
	@split_col8=split(";",$split_lines[7]);
	$isnotheader=substr($lines,0,1);
	@split_samplelist=split(",",$samplelist);	
	if ($isnotheader ne "#")
	{	
		# For each sample, comma-separated list of values
		for (my $e=0;$e<=$#split_samplelist;$e++)
		{
			for (my $i=0; $i<=$#split_col8; $i++)
			{
				$valCol=$split_col8[$i];
				$valCol=~s/=.*//;
				$value=$split_col8[$i];
				$value=~s/.*=//;
				if ($valCol =~ /$split_samplelist[$e]_D/)
        		{
					$depth=$value;
					if ($depth eq "-")
					{
						$depth=0;
					}
        		}
				elsif($valCol =~ /$split_samplelist[$e]_F/)
				{
					$allelefreq=$value;
				}
				elsif($valCol =~ /$split_samplelist[$e]_G/)
				{
					if (($value eq "P_Homo_ref") || ($value eq "Homo_ref"))
					{
						$genotype="0/0";
					}
					elsif (($value eq "UNC_Hetero") || ($value eq "P_Hetero"))
					{
						$genotype="0/1";
					}
					elsif (($value eq "P_Homo_var") || ($value eq "Homo_var"))
					{
						$genotype="1/1";
					}
					
					$ecRef=int($depth*(1-$allelefreq));
					$ecAlt=int($depth*$allelefreq);
					$depth=$ecRef+$ecAlt;
					push(@sampleInfo,$genotype.":".$ecRef.",".$ecAlt.":".$depth);
				}
			}
		}
		print VCF $split_lines[0],"\t",$split_lines[1],"\t",$split_lines[2],"\t",$split_lines[3],"\t",$split_lines[4],"\t",$split_lines[5],"\t",$split_lines[6],"\t.\tGT:EC:DP\t",join("\t",@sampleInfo),"\n";
		@sampleInfo=();
	}
}
close (PSEUDOVCF);
close (VCF);
