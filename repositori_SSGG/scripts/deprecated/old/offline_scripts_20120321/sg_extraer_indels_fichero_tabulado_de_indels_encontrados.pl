#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Extrae los indels identificados en la muestra que pertenecen a las zonas capturadas por el kit de Agilent (tres primeras columnas de un fichero bed) \n COMO SE USA:\n";
    print "Input: 1) fichero tabulado con las regiones que estudiamos. La primera columna es el cromosoma, la segunda la coordenada de inicio y la tercera, la coordenada final. Debe estar ordenado de menor a mayor.\nejemplo: chr1    343432      343445\n";
    print "       2) fichero tabulado con los indels en formato gff del bioscope\n";
    print " Output: fichero tabulado con los indels dentro de las coordenadas del exoma. Se hace por cromosoma y darle un fichero de salida. \n";
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
	#@split_lines2=split ("\t",$lines2);
	
	#$longitud_start=length($split_lines2[2]);
	#$longitud_final=length($split_lines2[3]);
	#if ($longitud_start > $longitud_final)
	#{
		push (@indels,$lines2);
	#}
	#elsif ($longitud_start <= $longitud_final)
	#{
	#	push (@indels,$lines2."\t".$longitud_final);
	#}
}
#chr10   deletion        245688  245693  6       245686  GTTCTTTATCTTGT/GTTCTTGT/NO_CALL REF,3,1 3       HOMOZYGOUS      0.0039
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
		$inicio_indel=$split_indels[2]-$split_indels[4];
		$final_indel=$split_indels[3]+$split_indels[4];
		#$inicio_indel=$inicio_indel-$split_indels[6];
		if (($inicio_indel >= $split_bedfile[1]) && ($inicio_indel <= $split_bedfile[2]) || (($inicio_indel <= $split_bedfile[1]) && ($final_indel >= $split_bedfile[1])) || (($inicio_indel <= $split_bedfile[2]) && ($final_indel >= $split_bedfile[2])))
		{
			print $indels[$l],"\n";
			$j++;
			$l=$l+1;
		}
		elsif($inicio_indel > $split_bedfile[2])
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
