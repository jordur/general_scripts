#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Template script for developing Perl scripts
# @Author: Bioinfo
# @Contributors: GOOD
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;




my $input=$ARGV[0];

my @explode=();
my @explode2=();
my @explode3=();

my %hash_samples=();

my $head_string="";
my $line;
my $string;
my $string2;
my $i;


open(FILE,"<",$input);
while($line=<FILE>)
{
	chomp($line);
	if($line=~/^#Samples: /)
	{
		$line=~s/#Samples: //;
		@explode=split(";",$line);
		for($i=0;$i<scalar(@explode);$i++)
		{
			@explode2=split("=",$explode[$i]);
			$hash_samples{$explode2[0]}=$explode2[1];
			$head_string=$head_string."\t".$explode2[1];
		}
		print "$head_string\n";
	}
	else
	{
		if($line=~/#/)
		{
			$head_string=$line;
		}
		else
		{
			$string="";

			@explode=split("\t",$line);
			for($i=0;$i<7;$i++)
			{
				if($string eq "")
				{
					$string=$explode[$i];
				}
				else
				{
					$string=$string."\t".$explode[$i];
				}
			}
			@explode2=split(";",$explode[7]);
			for($i=0;$i<scalar(@explode2);$i++)
			{
				@explode3=split("=",$explode2[$i]);
				if(!$explode3[1] || $explode3[1] eq "-")
				{
					$explode3[1]=0;
				}

				$string2=$explode3[0];
				if($explode3[0] =~/D/)
				{
					$string2=~s/D//;
					if($i==0)
					{
						$string=$string."\t".$hash_samples{$string2}."_D=".$explode3[1];
					}
					else
					{
						$string=$string.";".$hash_samples{$string2}."_D=".$explode3[1];
					}
				}
				elsif($explode3[0] =~/F/)
				{
					$string2=~s/F//;
					if($i==0)
					{
						$string=$string."\t".$hash_samples{$string2}."_F=".$explode3[1];
					}
					else
					{
						$string=$string.";".$hash_samples{$string2}."_F=".$explode3[1];
					}
				}
				elsif($explode3[0] =~/G/)
				{
					$string2=~s/G//;
					if($i==0)
					{
						$string=$string."\t".$hash_samples{$string2}."_G=".$explode3[1];
					}
					else
					{
						$string=$string.";".$hash_samples{$string2}."_G=".$explode3[1];
					}
				}
				else
				{
					$string=$string.";".$explode2[$i];
				}
			}
			$string=$string."\t".$explode[8]."\t".$explode[9];
			print "$string\n";
		}
	}
}
