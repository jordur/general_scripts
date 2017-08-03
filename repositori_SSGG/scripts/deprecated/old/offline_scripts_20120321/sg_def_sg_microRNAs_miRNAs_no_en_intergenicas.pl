#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: comprueba que ningun miRNA predicho por miRDeep está dentro de un gen\n";
    print " Input: 1) fichero parseado de la salida de miRDeep, 2) fichero .gff con los genes anotados, 3) nombre del cromosoma ej: 2L, 3R y 4) nombre del fichero output ej: CHR1 \n\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(MIRNAS,"<",@ARGV[0]);
open(GENES,"<",@ARGV[1]);
my $chr= @ARGV[2];
my $nombre_chr=@ARGV[3];
open(GENSENSE,">","$nombre_chr.gen.sense") or die "No puedo abrir $nombre_chr.gen.sense\n";
open(GENANTISENSE,">","$nombre_chr.gen.antisense") or die "No puedo abrir $nombre_chr.gen.antisense\n";
open(INTERGENICO,">","$nombre_chr.intergenico") or die "No puedo abrir $nombre_chr.intergenico\n";

my @mirnas;
my @genes;
my @split_linea;
my $contador=0;
while  (my $linea2=<MIRNAS>)
{
        chomp($linea2);
        #@split_linea= split ("_",$linea2);
	#$posicion_inicio=substr ($split_linea[2],2,length $split_linea[2]);
	#$signo= substr ($split_linea[2],0,1);
	#$nueva_linea = $linea2."\t".$posicion_inicio."\t".$split_linea[2]."\t".$signo;	
	#print $nueva_linea,"\n";
	push (@mirnas,$linea2);
}

close (MIRNAS);


while  (my $linea=<GENES>)
{
        chomp($linea);
        push (@genes, $linea);
}

close (GENES);


for (my $i =0; $i<= $#mirnas; $i++)
{

	$almohadilla= substr ($genes[$a],0,1);
	if ($almohadilla eq "#")
        {
	         $i++;
        }

	@split_linea= split ("\t",$mirnas[$i]);
		
	for (my $a=0; $a<=$#genes; $a++)
	{
		@split_linea_genes= split ("\t",$genes[$a]);
		@split_col9= split (";",$split_linea_genes[8]);
		#print $split_linea [7],"\t", $split_linea [8],"\t",$split_linea [9],"\n";
		if (($chr eq $split_linea_genes[0]) && ($split_linea [7] >= $split_linea_genes [3]) && ($split_linea [8] <= $split_linea_genes [4]) && ($split_linea [9] eq $split_linea_genes [6]))
		{
			print GENSENSE $mirnas[$i],"\t",$split_col9[0],"\t",$split_col9[1],"\n";
			$contador=$contador+1;
		}
		elsif (($chr eq $split_linea_genes[0]) && ($split_linea [7] >= $split_linea_genes [3]) && ($split_linea [8] <= $split_linea_genes [4]) && ($split_linea [9] ne $split_linea_genes [6]))
		{
			print GENANTISENSE $mirnas[$i],"\t",$split_col9[0],"\t",$split_col9[1],"\n";
			$contador=$contador+1;
		}
		
	}

	if ($contador == 0)
	{
		print INTERGENICO $mirnas[$i],"\n";
	}

	$contador=0;
}

close (SENSE);
close (ANTISENSE);
close (INTERGENICO);
