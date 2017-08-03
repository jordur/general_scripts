#!/usr/bin/perl -w

use POSIX qw/ strftime /;

open(WARR,'>',"errors.log");


if(!$ARGV[0])
{
	print WARR "Error\tPreprocess\tDon't find input\t";
	die;
}
$input=$ARGV[0];

my $count=0;
my $columns_comtrol=0;

open(FILE,'<',$input);
while($line=<FILE>)
{
	my $string="";
	my $scalar=0;

        chomp($line);
	if($line ne "")
	{
		if($line!~/^\/\//)
		{
			@explode=split("\t",$line);
			$scalar=scalar(@explode);
			if($count==0)
			{	
				$columns_comtrol=$scalar;
			}
			else
			{
				if($scalar != $columns_comtrol)
				{
					print WARR "Error\tPreprocess\t$count  Line has different columns's number\n";
				}
			}
			for(my $i=0;$i<scalar(@explode);$i++)
			{
				if($explode[$i] eq "")
				{
					$explode[$i]="-";
				}
				if($i==0)
				{
					$string="$explode[$i]";
				}
				else
				{
					$string=$string."\t".$explode[$i];
				}
			}
			print "$string\n";
			$count++;
		}
	}
}
