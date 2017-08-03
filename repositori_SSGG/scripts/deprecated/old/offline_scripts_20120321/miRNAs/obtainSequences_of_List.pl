#!/usr/bin/perl w-

$inputList=$ARGV[0];
$inputSequences=$ARGV[1];

%list=();

$entry=0;

open(FILE,'<',$inputList);
while($line=<FILE>)
{
	chomp($line);

	$list{$line}=1;
}
open(FILE2,'<',$inputSequences);
while($line=<FILE2>)
{
        chomp($line);

	if($line=~/^>/)
	{
		$lineGhost=$line;
		$lineGhost=~s/>//;
		if($list{$lineGhost}==1)
		{
			print "$line\n";
			$entry=1;
		}
		else
		{
			$entry=0;
		}
	}
	else
	{
		if($entry==1)
		{
			print "$line\n";
			$enrty=0;
		}
	}
}
