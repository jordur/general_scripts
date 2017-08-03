#!/usr/bin/env perl
########################################################
# Copyright 2014 - Sistemas Genomicos S.L.
# @Desc: Class definition and some functions to map data
# @Author: Sheila
# @Contributors: 
########################################################

# -----------------------
# -- Include libraries --
# -----------------------

use warnings;
use strict;


# -----------------------
# -- Global variables  --
# -----------------------


my $linepileup ="chr17	41249363	T	8017	...-1A...-1A...-1A.";

&main;

sub count_reads_in_strand
{
    my ($linepileup) = @_;
    #print $linepileup,"\n";
    my @splitlinepileup=split("\t",$linepileup);
    my @splitreadstrands;
    my $positivestrand=0;
    my $negativestrand=0;
    my $refallele=0;
    my $alternativeallele=0;
    
    for (my $i=4;$i<=$#splitlinepileup;$i=$i+3)
    {
        print  $splitlinepileup[0],"\t",$splitlinepileup[1],"\t";
        @splitreadstrands=split("",$splitlinepileup[$i]);
        for (my $a=0;$a<=$#splitreadstrands;$a++)
        {
        	#print $splitreadstrands[$a],"\n";
        	if (($splitreadstrands[$a] =~ /\./) || ($splitreadstrands[$a] =~ />/))
        	{
        		$positivestrand=$positivestrand+1;
        	}
        }
        print $positivestrand,"\n";

    }
}

sub main
{
	&count_reads_in_strand($linepileup);
}
