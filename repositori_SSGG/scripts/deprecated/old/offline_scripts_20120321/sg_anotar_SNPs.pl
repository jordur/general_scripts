#!/usr/bin/perl

open(FILE1,"<",@ARGV[0]);
open(FILE2,"<",@ARGV[1]);
my @fichero1;
my @fichero2;
my @fila1;
my @fila2;
for (my $i=1; my $lines=<FILE1>; $i++)
{
  chomp($lines);
  push (@fichero1,$lines);
}

close (FILE1);


for (my $i=1; my $lines2=<FILE2>; $i++)
{
  chomp($lines2);
  push (@fichero2,$lines2);
}

close (FILE2);


for (my $i=0; $i<=$#fichero2; $i++)
{
	@fila2 = split (/\t/,$fichero2[$i]);
	for (my $i=0; $i<=$#fichero1; $i++)
	{
		
		@fila1 = split (/\t/,$fichero1[$i]);
		if ((@fila2[11] >= $fila1[3]) && (@fila2[11] <= $fila1[4]))
		{
			@Gene_ID=split (";",$fila1[8]);
			#$posicion = index ("protein_id=",4);
			#print $posicion,"\n"; 
			#print @Gene_ID;	
			print $fila2[0],"\t",$fila2[1],"\t",$fila2[2],"\t",$fila2[3],"\t",$fila2[4],"\t",$fila2[5],"\t",$fila2[6],"\t",$fila2[7],"\t",$fila2[8],"\t",$fila2[9],"\t",$fila2[10],"\t",$fila2[11],"\t",$fila1[3],"\t",$fila1[4],"\t",$Gene_ID[0],"\t",$fila1[8],"\n";
		}
		#elsif ((@fila2[11] != $fila1[3]) && (@fila2[11] !> $fila1[3]) && (@fila2[11] != $fila1[4]) && (@fila2[11] !< $fila1[4]))	
		#{
		#	print $fila2[0],"\t",$fila2[1],"\t",$fila2[2],"\t",$fila2[3],"\t",$fila2[4],"\t",$fila2[5],"\t",$fila2[6],"\t",$fila2[7],"\t",$fila2[8],"\t",$fila2[9],"\t",$fila2[10],"\t",$fila2[11],"\n";		
		#}
	}
}

