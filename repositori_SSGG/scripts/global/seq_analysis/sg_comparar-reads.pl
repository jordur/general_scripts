#!/usr/bin/perl -w

######################

## Date: 3/8/12

######################



use Bio::SeqIO;



sub usage
{
        print " USO: Parte del pipe sg_colapsar-lecturas-paired-end.sh. Compara las lecturas que han sido previamente ordenadas en un paso anterior para generar una unica lectura no redundante en formato fasta en cuya cabecera se puede ver el numero de veces que se repite en el dataset \n";
        print " Input: ";
        print "         \$1: fichero temporal generado por sg_collapse-reads.pl \n";
        print " Output: Dos ficheros, uno para cada grupo de lecturas, con aquellas lecturas no redundantes \n";
        exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}



open (FILE,"<",$ARGV[0]);
open (READ1,">","p1-read1.fastq");
open (READ2,">","p1-read2.fastq");
$totalLineas=$ARGV[1];
$contador=0;
$exp=1;
$lecturas=0;

while ($lines=<FILE>)
{
	chomp ($lines);
	$contador=$contador+1;
	if ($contador > 1)
	{
		@split_linea=split("\t",$linea);
		if ($contador < $totalLineas && $lines eq $linea)
		{
			$exp=$exp+1;
		}
		elsif ($contador < $totalLineas && $lines ne $linea)
		{
			$lecturas=$lecturas+1;
			print READ1 ">",$lecturas."_".$exp,"\n",$split_linea[0],"\n";
			print READ2 ">",$lecturas."_".$exp,"\n",$split_linea[1],"\n";
			$exp=1;
		}
		elsif ($contador == $totalLineas && $lines ne $linea)
		{
			$lecturas=$lecturas+1;
			print READ1 ">",$lecturas."_".$exp,"\n",$split_linea[0],"\n";
			print READ2 ">",$lecturas."_".$exp,"\n",$split_linea[1],"\n";
			$lecturas=$lecturas+1;
			print READ1 ">",$lecturas."_1\n",$split_linea[0],"\n";
			print READ2 ">",$lecturas."_1\n",$split_linea[1],"\n";
		}
		elsif ($contador == $totalLineas && $lines eq $linea)
                {
                        $exp=$exp+1;
			$lecturas=$lecturas+1;
                        print READ1 ">",$lecturas."_".$exp,"\n",$split_linea[0],"\n";
			print READ2 ">",$lecturas."_".$exp,"\n",$split_linea[1],"\n";
                }

	}
	$linea=$lines;
	
}
