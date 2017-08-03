#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los indels de un fichero tabulado que se encuentran en las regiones indicadas por un archivo donde indicamos los intervalos de secuencias que estamos analizando, tres primeras columnas de un fichero bed \n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) fichero tabulado con los indels en el que la primera columna es el cromosoma\n";
    print " Output: fichero tabulado con los SNPs que se encuentran en las regiones indicadas en el fichero de input 1. Hay que hacerlo cromosoma por cromosoma y darle un fichero de salida. \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(BEDFILE,"<",@ARGV[0]);
open(INDELS,"<",@ARGV[1]);
my @bedfile;
my @indels;
my @split_bedfile;
my @split_indels;

while ($lines2=<INDELS>)
{
	chomp($lines2);
	@split_lines2=split ("\t",$lines2);
	
	$longitud_start=length($split_lines2[2]);
	$longitud_final=length($split_lines2[3]);
	if ($longitud_start > $longitud_final)
	{
		push (@indels,$lines2."\t".$longitud_start);
	}
	elsif ($longitud_start <= $longitud_final)
	{
		push (@indels,$lines2."\t".$longitud_final);
	}
}

close (INDELS);

my $l = 0;
while (my $bed=<BEDFILE>)
{
	chomp ($bed);
	

        @split_bedfile = split (/\t/,$bed);
	for (my $j=0; $j<=$#indels;)
	{
		$j=$l;
		@split_indels = split (/\t/,$indels[$l]);
		$final_indel=$split_indels[1]+$split_indels[6]+$split_indels[6];
		$split_indels[1]=$split_indels[1]-$split_indels[6];
		if (($split_indels[1] >= $split_bedfile[1]) && ($final_indel <= $split_bedfile[2])|| (($split_indels[1] <= $split_bedfile[1]) && ($final_indel >= $split_bedfile[1])) || (($split_indels[1] <= $split_bedfile[2]) && ($final_indel >= $split_bedfile[2])))
		{
			print $indels[$l],"\n";
			$j++;
			$l=$l+1;
		}
		elsif($split_indels[1] > $split_bedfile[2])
		{	
			
			last;
			
		}
		else 
		{
			$j++;
			$l=$l+1;
		}
		
	}
}
close(BEDFILE);
