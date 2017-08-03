#!/usr/bin/env perl
########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map Illumina data using bwa
# @Author: JM Rosa, Arbol
# @Contributors: 
########################################################

package mapping_bioscope;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use NGS::essay;
use file_handle;
use sam_handle;
use bam_handle;
use submit;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------
my $module = "mapping";
our @ISA = qw(Exporter);
our @EXPORT = qw();

sub mapping{
	my ($class, $essay, $jobs) = @_;
	print "Launching mapping_bioscope jobs for essay $essay->{name}... ".localtime()."\n";
	my @types = ("F3", "PE");
	my $list3;
	foreach my $type (@types){
		my $list2;
		foreach my $sample_name (keys(%{$essay->get_parameter("samples")})){
			foreach my $replicate_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates"))})){
				my $list;
				foreach my $lane_name (keys(%{$essay->get_parameter(("samples",$sample_name,"replicates",$replicate_name,"lanes"))})){
					my $bam_id = sam_handle->manage_sam($essay,$jobs,$sample_name,$replicate_name,$lane_name,$type);
					$list = join_jobs($list,$bam_id);
				}
				my @files = @{$essay->get_bams(sample=>$sample_name,replicate=>$replicate_name,level=>"lane",type=>$type)};
				my $merge_id = &merge_lanes (essay=>$essay,jobs=>$list,files=>\@files,sample=>$sample_name,replicate=>$replicate_name,type=>$type);
				my $hq_bam_id = &get_hq_bam ($essay,$merge_id,$sample_name,$replicate_name,$type);
				my $realign_id = &realign ($essay,$hq_bam_id,$sample_name,$replicate_name,$type);
				$list2 = join_jobs($list2,$realign_id);
			}
		}
		#my @files = @{$essay->get_bams(type=>$type)};
		#my $merge_id = &merge_lanes (essay=>$essay,jobs=>$list2,files=>\@files,type=>$type);
		my $list3 = join_jobs($list3,$list2);
	}
	return $list3;
}

1;