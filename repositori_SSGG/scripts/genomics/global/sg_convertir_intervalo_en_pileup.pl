#!/usr/bin/perl -w
#use strict;

my @split_linea;

open (INTERVALOS,"<",$ARGV[0]);
while (my $lineas = <INTERVALOS>)
{
	chomp ($lineas);
	@split_linea = split ("\t",$lineas);
	for (my $i=$split_linea[1]; $i<=$split_linea[2]; $i++)
	{
		print "$split_linea[0]\t$i\t-\t-\t-\t-\t-\t20\n";
	}
	
}
close (INTERVALOS);
