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
#use MIME::Base64;

# Debugging (optional)
#SOAP::Lite->import(+trace => qw(debug));

#my @genes;

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

#open GENES, "3gene.txt";
#
#while (<GENES>){
#	chomp $_;
#	push(@genes, $_);	
#}

#	foreach my $gene (@genes){

my $argsXrefData = SOAP::Data->name(
	ids => 'ENSG00000139618',
	codes => 'En',
	);
	
my @foundPathwayRefList = $wp_soap->findPathwaysByXref($argsXrefData)->paramsout;
unshift(@foundPathwayRefList, $wp_soap->findPathwaysByXref($argsXrefData)->result); #add first result to list
	#print "\nfindPathwaysByXref(id=$gene, code=X):\n";
	foreach my $ref (@foundPathwayRefList){
		#print $ref,"\n";
		print "$ref->{id}\n";
		foreach my $key (keys %{ $ref}){
			if ($key =~ /^fields$/){ # Print all fields
				my $argsgetPathwayData = SOAP::Data->name(
					pwId => $ref->{id},
					revision => 0,
				);
				
				my $pw_res = $wp_soap->getPathway($argsgetPathwayData)->result;
				
				foreach my $f (@{ $ref->{$key} }) {
					my $name = $f->{'name'};
					my $value = $f->{'values'};
					print "$name: $value\n";
				}
			} else {
				#print "\t$key: $ref->{$key}\n";
			}
		}
		print "############\n\n";
	}
