#!/usr/bin/env perl
##############################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to launch Miseq Pipeline with BWA
# @Desc: Miseq enhanced Launcher in Simple Appareance
# @Author: JM Rosa
# @Contributors: Arbol
###############################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use file_handle;


# ---------------------------
# -- functions definitions --
# ---------------------------


# ---------------------------------
# -- global & default parameters --
# ---------------------------------
my $metadata;
my $cwd = `pwd`;
chomp $cwd;
my $output_id = cwd() . "/output_file";

# --------------------------------
# -- Input parameters & options --
# --------------------------------
(my $basename, my $dir, my $ext) = fileparse($0, qr/\.[^.]*/);

GetOptions (
				"m=s"	=> \$metadata,
				"o=s"	=> \$output_id,
);

print "\n=================================================================
$basename: Miseq enhanced Launcher in Simple Appareance v1.0\n";

if (not defined $metadata) {
	die "
			Options:
		-m Metadata File
		-o Ouput name (optional)
\n".localtime()."\n=================================================================\n\n";
}

print  "\nProcess started at... ".localtime()."\n";
&main ();
print  "Process finished... ".localtime()."\n";

# ---------------
# -- functions --
# ---------------
sub main (){

	print "Running $basename with parameters m=$metadata and o=$output_id\n";
	
	&launch_json2essay ($metadata, $output_id);
	
	my $essays_file = &get_in_file_handle ($output_id.".tsv");
	
	while (my $line = <$essays_file>) {
		chomp $line;
		if ($line =~ m/^#/) {next;}
		my $essay = (split(/\s+/, $line))[0];
		my $json = (split(/\s+/, $line))[1];
		my $aligner = (split(/\s+/, $line))[2];
		my $mapping_id = &get_mapping_id ($essay, $json, $aligner); 
		my $var_calling_id = &launch_variant_calling  ($essay, $json, $mapping_id);
		my $anno_id = &launch_annotation  ($essay, $json, $var_calling_id);
		print "$anno_id\n";
	}
}

sub get_mapping_id {
	my ($self, $file, $mode) = @_;
	if ($mode eq 'bwa'){
		return &launch_mapping_bwa($self, $file);
	}	
	elsif ($mode eq 'novoalign'){
		return &launch_mapping_novoalign($self, $file);
	}
	else { die "Aligner not recognised\n";}
}

sub launch_annotation {
	my ($self, $file, $id) = @_;
	print "Launching annotation.pl -m $file -o $self -hold_jid $id > annotation_$self.log... ".localtime()."\n";
	`annotation.pl -m $file -o $self -hold_jid $id > annotation_$self.log`;
	#&check_log (variant_$self.log);
	my $col = '$2';
	my $job = `grep JOB_ID annotation_$self.log | awk '{print $col}'`;
	chomp $job;
	return $job;		
}

sub launch_variant_calling {
	my ($self, $file, $id) = @_;
	print "Launching variant_calling.pl -m $file -o $self -hold_jid $id > variant_$self.log... ".localtime()."\n";
	`variant_calling.pl -m $file -o $self -hold_jid $id > variant_$self.log`;
	#&check_log (variant_$self.log);
	my $col = '$2';
	my $job = `grep JOB_ID variant_$self.log | awk '{print $col}'`;
	chomp $job;
	return $job;		
}

sub launch_mapping_novoalign {
	my ($self, $file) = @_;
	print "Launching mapping_novoalign.pl -m $file -o $self > mapping_$self.log... ".localtime()."\n";
	`mapping_novoalign.pl -m $file -o $self > mapping_$self.log`;
	#&check_log (mapping_$self.log);
	my $col = '$2';
	my $job = `grep JOB_ID mapping_$self.log | awk '{print $col}'`;
	chomp $job;
	return $job;		
}

sub launch_mapping_bwa {
	my ($self, $file) = @_;
	print "Launching mapping_bwa.pl -m $file -o $self > mapping_$self.log... ".localtime()."\n";
	`mapping_bwa.pl -m $file -o $self > mapping_$self.log`;
	#&check_log (mapping_$self.log);
	my $col = '$2';
	my $job = `grep JOB_ID mapping_$self.log | awk '{print $col}'`;
	chomp $job;
	return $job;		
}

sub launch_json2essay {
	my ($json, $output) = @_;
	print "json2essay.pl -m $json -o $output > $json.log... ".localtime()."\n";
	`json2essay.pl -m $json -o $output > $json.log`;
	#&check_log ($json.log);
	return 1;
}

exit;