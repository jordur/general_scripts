#!/usr/bin/perl -w

#El input es el fichero de los miRNAs maduros de mirBase. El objetivo de este script es separar las secuencias de miRNAs maduros de las miRNA star

open (IN,"<",$ARGV[0]);
open (MADURO,">","miRNAs_maduros.fasta");
open (STAR,">","miRNAs_star.fasta");
$contador=0;
while ($linea = <IN>)
{
	chomp ($linea);
	$primero = substr ($linea,0,1);
	if (($primero  eq ">") && (asterisco($linea) eq "No"))
	{
		$contador=1;
		print MADURO $linea,"\n";
	}
	elsif (($primero ne ">") && ($contador == 1))
	{
		print MADURO $linea,"\n";
		$contador=0;
	}
	else
	{
		print STAR $linea,"\n";
	}
	
}

close (IN);
close (MADURO);
close (STAR);

sub asterisco
{

	$string = shift;
	if($string =~m/\*/)
	{ 
		return "Yes"; 
	}
	else
	{
		return "No";
	}
}
