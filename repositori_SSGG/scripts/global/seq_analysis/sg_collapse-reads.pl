#!/usr/bin/perl -w

######################

## Date: 3/8/12


######################



use Bio::SeqIO;



sub usage
{
        print " USO: Parte del pipe sg_colapsar-lecturas-paired-end.sh. Formatea un fichero temporal que contiene el resultado de pegar los dos archivos que contienen las lecturas con los paired-ends  \n";
        print " Input: ";
        print "         \$1: un unico fichero con las lecturas paired end en la misma linea \n";
        print " Output: Fichero temporal formateado \n";
        exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}



open (READ1,"<",$ARGV[0]);
$contador=0;

while ($linea1=<READ1>)
{
	chomp ($linea1);
	$contador=$contador+1;
	if ($contador == 2)
	{
		print $linea1,"\n";
		
	}

	elsif ($contador==4)
	{
		$contador=0;		
	}
}
