#! usr/bin/perl -w
use strict;

########## Script per seleccionar/desseleccionar un sèrie de valors d'un array d'una altra sèrie de valors.(en aquest cas desseleccionar les id de serps del resultat blast)###################
open FILE1, "<all_hits.txt";
open FILE2, "<all_gi_serpentes";
open FILE3, ">output";
my @hit = <FILE1>;
my @serp= <FILE2>;
my %id;
# print "@hit\n";
# print "@serp\n";
foreach my $line ( @serp ){
 chomp $line;
 $id {$line} = 1;	
}

# print keys%id;

my @collect_lines;

foreach my $b (@hit){
# print "$b";	
	if ( my ($gi) = $b =~ /^.+gi\|(\d+)/) {
		if (! $id { $gi } ){
			push ( @collect_lines, $b );
# 		} else {
# 			print " I don't want this gi $gi\n";
		}
# print @a;
	} else {
		print "this file has no gi\n";
	}
}

print FILE3 "@collect_lines";



############## Una forma similar de fer-ho és usant map (no ho he provat)#################
# %id = map {$_, 1} @serp;
# foreach my $b (@hit){
# if ( my ($gi) = $b =~ /^.+gi\|(d+)/) {
# my @difference = grep {!$id {$gi}} $b;
# 
# print "@difference\n";
# 	}
# }	
