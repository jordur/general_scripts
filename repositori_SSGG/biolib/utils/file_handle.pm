#!/usr/bin/env perl
###############################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to work with files
# @Author: JM Rosa
# @Contributors:
###############################################################

package file_handle;

# -----------------------
# -- Include libraries --
# -----------------------
use warnings;
use strict;
use FileHandle;
use Exporter;

# -----------------------
# -- Global variables  --
# -----------------------

our @ISA = qw(Exporter);
our @EXPORT = qw (get_in_file_handle get_out_file_handle);

sub get_in_file_handle {
	
	my $file = shift;
	
    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($file)) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $file, "\n") unless -e $file;
        
        $in_file_handle->open( $file ) or die("ERROR: Could not read from input file ", $file, "\n");
    }
    
    # no file specified - try to read data off command line
    else {
		print "INFO: Reading input from STDIN (or maybe you forgot to specify an input file?)...";
		$in_file_handle = 'STDIN';
    }
    
    return $in_file_handle;
}

sub get_out_file_handle {
	my ($file, $extension) = @_;
	
	if (defined($extension)){
		$extension = ".".$extension;
	} else {
		$extension = "";
	}
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open("> ".$file.$extension) or die("ERROR: Could not write to output file ".$file.$extension."\n");

	return $out_file_handle;
}

1;