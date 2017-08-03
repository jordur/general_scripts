#!/usr/bin/env perl
########################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Library containing transcripts utilities
# @Author: Arbol
# @Contributors: biomart & ensembl scripts
########################################################

package transcripts;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(new_transcript show_transcript);

# -----------------------
# -- Include libraries --
# -----------------------

# ---------------
# -- functions --
# ---------------

sub new_transcript {
	# This is an auxiliar subroutine for creating new transcripts
	
	# Arguments:
	my ( $transcript, $strand ) = @_;
	
	${$transcript}{'size'} = 0;
	${$transcript}{'UTR5'} = 0;
	${$transcript}{'UTR3'} = 0;
	${$transcript}{'strand'} = $strand;
	${$transcript}{'found_UTR5'} = 0;
	${$transcript}{'found_UTR3'} = 0;
	${$transcript}{'exons'} = ();
	${$transcript}{'UTR'} = ();
}

sub show_transcript {
	# This is an auxiliar subroutine for showing transcripts
	
	# Arguments:
	my ( $transcript ) = @_;
	
	foreach ( keys %{$transcript} ) {
		print $_ . ": " . ${$transcript}{$_} . "\n";
	}
}


# **************************************************************************************************
# **************************************************************************************************
1;
 