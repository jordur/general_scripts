#!/usr/bin/perl -w

use strict;

if($ARGV[0] eq "")
{
	print "\n\n\t\tI:Input File (promer file with fragment Genome\n";
}
else
{

	my $header="";
	my $input=$ARGV[0];
	my $line="";
	
	open(FILE,'<',$input);
	while($line=<FILE>)
	{
		chomp($line);
		if($line=~(/>/))
		{
			$header=$line;
		}
		else
		{
			$header=~s/>//gi;




			print "$header\t$line\n";
		}
	}
}
