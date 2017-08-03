#!/usr/bin/perl

my $input=$ARGV[0];

my @explode=();

my $line="";

open(FILE,"<",$ARGV[0]);
while($line=<FILE>)
{
	chomp($line);

	if($line=~/\w/)
	{
		@explode=split("\t",$line);
		for(my $i=$explode[1];$i<($explode[2]+1);$i++)
		{
			print "$explode[0]\t$i\t.\n";
		}
	}
}
