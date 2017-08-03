#!/usr/bin/perl -w
use Math::NumberCruncher;

#my $mode = $ARGV[0];
open (COVERAGE,"<",$ARGV[0]);
	while (my $lineas = <COVERAGE>)
	{
		chomp ($lineas);
		#if ($mode eq "columna")
		#{
		@divido_linea=split ("\t",$lineas);
		for (my $i=2;$i<=$#divido_linea;$i++)
		{
			push (@numArray,$divido_linea[$i]);
			
		}
		$median = Math::NumberCruncher::Median(\@numArray,1);
		print $divido_linea[0],"\t",$divido_linea[1],"\t";
		for (my $a=2;$a<$#divido_linea;$a++)
		{
			if ($divido_linea[$a] == 0 || $median == 0)
			{
				print "0\t";
			}
			else
			{
				print $divido_linea[$a]/$median,"\t";
			}
		}
		if ($divido_linea[$#divido_linea] == 0 || $median == 0)
		{
			print "0\n";
		}
		else
		{
			print $divido_linea[$#divido_linea]/$median,"\n";
		}
		@numArray=();
	}


close (COVERAGE);
