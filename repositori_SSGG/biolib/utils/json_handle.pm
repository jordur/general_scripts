#!/usr/bin/env perl
###############################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with files
# @Author: JM Rosa & Arbol
# @Contributors:
###############################################################

package json_handle;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use JSON -support_by_pp;
use JSON -convert_blessed_universally;
use file_handle;
use Exporter;
use Data::Dumper;

# -----------------------
# -- Global variables  --
# -----------------------

our @ISA = qw(Exporter);
our @EXPORT = qw (get_from_json trasform_to_json create_json);

sub get_from_json {
	my ($file) = @_;
	
	my $json;
	{
		local $/; #Enable 'slurp' mode
		open my $fh, "<", $file or die "ERROR: Unable to open JSON file $file!!\n";
		$json = <$fh>;
		close $fh;
	}
	my %data = %{decode_json($json)};
    
    return %data;
}

sub trasform_to_json {
	my $self = shift;
	return JSON->new->canonical->allow_blessed->convert_blessed->encode($self);  
}

sub create_json {
	my ($object,$file) = @_;
	my $json_string = JSON->new->pretty->allow_blessed->convert_blessed->encode($object);
	my $output = get_out_file_handle($file, 'json');
	print $output $json_string;
	return $file.".json";
}

1;