#!/usr/bin/perl -w

open (FILE,"<",$ARGV[0]);
while (my $lines=<FILE>)
{
	chomp($lines);
	$char= substr($lines,0,1);
	if ($char eq ">")
	{
		$lines=~s/>//;
		$contigName=$lines;
	}
	else
	{
		@split_line=split(" ",$lines);
		print $split_line[0],"\t",$contigName,"\t",$split_line[1],"\t",$split_line[2],"\t",$split_line[3],"\n";
	}
}
