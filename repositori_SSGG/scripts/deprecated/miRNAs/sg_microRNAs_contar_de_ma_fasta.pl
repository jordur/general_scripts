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

#print join ("\n",@tab),"\n";
my $contador=1;
my $a;

for (my $i=0; $i<=$#tab;$i++)
{
	$a=$i+1;
        @split_tab1 = split (/\t/,$tab[$i]);
        @split_tab2 = split (/\t/,$tab[$a]);
        $valor_inferior = $split_tab1[1];
        $valor_problema = $split_tab2[1];
	$valor_superior = $split_tab1[2];
	#print $valor_superior,"\n";
	$valor_problema_superior = $split_tab2[2];
	#print $valor_superior,"\n";
	#print $valor_problema_superior,"\n";
	#print $valor_inferior,"\n";
	#push (@maximo, $valor_superior);
	#print $split_tab1[3],"\t",$split_tab2[3],"\n";	
	if (($a<=$#tab) && ((($valor_problema-$valor_inferior)<5)  || ($valor_inferior == $valor_problema)))
        {
		#print $split_tab1[3],"\t",$split_tab2[3],"\n";
		push (@number,$valor_inferior);
		push (@number,$valor_problema);
		push (@number,$valor_superior);
		push (@number,$valor_problema_superior);
		
		$contador=$contador+1;
	}
	else
	{
		if ($contador >=10)
		{
			$max = $number[0];
    			$min = $number[0];

    			foreach $i (@number[1..$#number])
    			{
        			if ($i > $max)
        			{
            				$max = $i;
        			}
        			elsif ($i < $min)
        			{
            				$min = $i;
        			}
    			}
			print $split_tab1[0]."_x".$contador,"\t",$min,"\t",$max,"\t",$split_tab1[3],"\n";
		}
		$contador=1;
		@number=();
	}
}
