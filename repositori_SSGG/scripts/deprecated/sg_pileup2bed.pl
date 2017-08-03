#!/usr/bin/perl

$input=$ARGV[0];##PILEUP INPUT BY CHR
@explode=();

$start=0;
$end=0;
$count=0;

open(FILE,"<",$ARGV[0]);
while($line=<FILE>)
{
        chomp($line);
	@explode=split("\t",$line);
	if($count==0)
	{
		$start=$explode[1];
		$end=$explode[1];
		$count=1;
	}
	else
	{
		if($explode[1]==$end+1)
		{
			$end=$explode[1];
		}
		else
		{
			print "$explode[0]\t$start\t$end\n";
			$start=$explode[1];
			$end=$explode[1];
		}
	}
}
if($start ne "")
{
	print "$explode[0]\t$start\t$end\n";
}
