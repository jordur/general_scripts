#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: Selecciona las sondas que caen dentro de los intervalos de genes que queremos analizar \n COMO SE USA:\n";
    print " Input: 1) fichero tabulado con 4 columnas, cromosoma, inicio gen, final del gen y nombre del gen \nejemplo: chr14   23881947        23904927        MYH7\n";
    print "        2) fichero tabulado de las sondas con 3 columnas: cromosoma, inicio de la sonda y final de la sonda \nejemplo: chr14    20138   20258\n";
    print " Output: sondas que se encuentran en los genes con el mismo formato del fichero tabulado de las sondas al que se le ha añadido una columna con el nombre del gen. Hay que dar un fichero de salida. Se debe hacer cromosoma a cromosoma. Si no coinciden las dos primeras columnas de ambos fichero, el programa da error. \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(GENES,"<",@ARGV[0]);
open(SONDAS,"<",@ARGV[1]);
my @genes;
my @sondas;
for (my $i=1; my $lines=<GENES>; $i++)
{
  chomp($lines);
  push (@genes,$lines);
}

close (GENES);

for (my $i=1; my $lines2=<SONDAS>; $i++)
{
  chomp($lines2);
  push (@sondas,$lines2);
}

close (SONDAS);

for (my $i=0; $i<=$#genes; $i++)
{
        @split_genes = split (/\t/,$genes[$i]);
	
        for (my $a=0; $a<=$#sondas; $a++)
        {
                @split_sondas = split (/\t/,$sondas[$a]);
                if($split_genes[0] eq $split_sondas[0])
		{
			if((($split_sondas[1] >= $split_genes[1]) && ($split_sondas[1] <= $split_genes[2])) || (($split_sondas[2] >= $split_genes[1]) && ($split_sondas[2] <= $split_genes[2])) )
			{
                        
				print $split_sondas[0],"\t",$split_sondas[1],"\t",$split_sondas[2],"\t",$split_genes[3],"\n";
                	}
        	}
		else
		{
			print "no coincide el cromosoma";
			last;
		}
	}
}

