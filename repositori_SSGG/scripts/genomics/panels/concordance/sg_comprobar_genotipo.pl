#!/usr/bin/perl -w

open (FILE,"<",$ARGV[0]);

while ($lineas=<FILE>)
{
	$alelomuestra="";
	$genotipo="";
	chomp($lineas);
	@splitLinea=split("\t",$lineas);
	@splitchar=split(//,$splitLinea[8]);
	if (($splitchar[0] eq $splitchar[1]) && ($splitchar[0] ne "N"))
	{
		$alelomuestra = $splitchar[0];
		if ($splitchar[0] eq $splitLinea[9])
		{
			$genotipo = "0/0"; 
		}
		elsif ($splitchar[0] ne $splitLinea[9])
		{		
			$genotipo = "1/1";
		}
	} 
	elsif (($splitchar[0] ne $splitchar[1]) && ($splitchar[0] ne "N"))
	{
		if ($splitchar[0] eq $splitLinea[9]) 
		{
			$alelomuestra =	$splitchar[1];
			$genotipo = "0/1";
		}
		elsif ($splitchar[1] eq $splitLinea[9])
		{
			$alelomuestra = $splitchar[0];
			$genotipo = "0/1";
		}
		elsif (($splitchar[0] ne $splitLinea[9]) && ($splitchar[1] ne $splitLinea[9]))
		{
			$alelomuestra = $splitchar[0].",".$splitchar[1];
			$genotipo = "1/2";
		}
	}
	
	if ($alelomuestra ne "")
	{
		print "$splitLinea[0]\t$splitLinea[3]\t$splitLinea[1]\t$splitLinea[9]\t$alelomuestra\t.\t.\t.\tGT:DP:EC\t$genotipo:.:.\n";
	}
	
}
