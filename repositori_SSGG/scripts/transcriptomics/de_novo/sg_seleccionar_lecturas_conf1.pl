#!/usr/bin/perl -w

sub usage
{
	print " USO: filtra lecturas de un fichero fasta atendiendo a su valor de confianza \n";
	print " Input: \n";
	print " 	\$1: Fichero fasta \n";
	print " 	\$2: Valor de confianza \n";
	
	print "	Output: Fichero fasta con las lecturas filtradas por su valor de confianza \n";
	exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


open (FICHERO,"<",$ARGV[0]);

while(my $lineas= <FICHERO>)
{
	chomp ($lineas);
	$char=substr ($lineas,0,1);
	if ($char eq ">")
	{
		if ($lineas =~ /Confidence_1.000/)
		{
			print $lineas,"\n";
			$imprimo="yes";
		}
		else
		{
			$imprimo="no";
		}
	}
		
	
	elsif ($imprimo eq "yes")
	{
		print $lineas,"\n";
		
	}

}

