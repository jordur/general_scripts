#!/usr/bin/perl
sub usage {
    print "\nPARA QUE SIRVE ESTO: Hace el conteo de las lecturas a patir del fichero .ma.tab\n";
    print "Input: fichero .ma.tab ordenado por posicion de inicio \n";
    print "Output: fichero en fasta con el conteo de las secuencias \n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(TAB,"<",@ARGV[0]);
my @tab;
my $contador=0;
while  (my $tab_linea=<TAB>)
{
  	chomp($tab_linea);
	push (@tab, $tab_linea);
}

close (TAB);
my $contador=1;
for (my $i=0; $i<=$#tab;)
{
	$a=$i+1;
        @split_tab1 = split (/\t/,$tab[$i]);
        @split_tab2 = split (/\t/,$tab[$a]);
        $valor_inferior = $split_tab1[1];
        $valor_problema = $split_tab2[2];
        if ((($split_tab2[1]-$valor_inferior)<5)  || ($valor_inferior == $split_tab2[1]))
        {
                $contador=$contador+1;
                $i++;
	}
	else
        {
                if ($contador >=2)
		{
			print $split_tab1[0]."_x".$contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		}
		$contador=1;
                $i++;
	}
}
