#!/usr/bin/perl -w
use Math::NumberCruncher;

#my $mode = $ARGV[0];
open (COVERAGE,"<",$ARGV[0]);
	while (my $lineas = <COVERAGE>)
	{
		chomp ($lineas);
		#if ($mode eq "columna")
		#{
		
			push (@numArray,$lineas);
		#}
		#elsif ($mode eq "linea")
	}


close (COVERAGE);

$median = Math::NumberCruncher::Median(\@numArray,1);
print $median,"\n";
