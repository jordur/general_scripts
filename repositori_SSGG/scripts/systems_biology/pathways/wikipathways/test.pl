#!/usr/bin/env perl

########################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: 
# @Author: Guille
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;

# Load module
use SOAP::Lite;

# Debugging (optional)
#SOAP::Lite->import(+trace => qw(debug));

# Create service interface with fault handler
my $wp_soap = SOAP::Lite
	->proxy('http://www.wikipathways.org/wpi/webservice/webservice.php')
	->uri('http://www.wikipathways.org/webservice')
	->on_fault(sub {
		my $soap = shift;
		my $res = shift;

		# Map faults to exceptions
		if(ref($res) eq '') {
			die($res);
		} else {
			die($res->faultstring);
		}
		return new SOAP::SOM;
} );

# findPathwaysByXref (multiple args; returns array of hash references)
my $argsXrefData = SOAP::Data->name(
	ids => '261_s_at',
	codes => 'X',
);

my @foundPathwayRefList = $wp_soap->findPathwaysByXref($argsXrefData)->paramsout;
unshift(@foundPathwayRefList, $wp_soap->findPathwaysByXref($argsXrefData)->result); #add first result to list
print "\nfindPathwaysByXref(id=205820_s_at, code=X):\n";
foreach my $ref (@foundPathwayRefList){
	foreach my $key (keys %{ $ref}){
		if ($key =~ /^fields$/){ # Print all fields
			foreach my $f (@{ $ref->{$key} }) {
				my $name = $f->{'name'};
				my $value = $f->{'values'};
				print "$name: $value";
			}
		} else {
			print "\t$key: $ref->{$key}\n";
		}
	}
}