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
        $valor_problema = $split_tab2[1];
	#print $split_tab1[0],"\t";
        if ((($valor_problema-$valor_inferior)<5)  || ($valor_inferior == $valor_problema))
        {
                $contador=$contador+1;
		#print $split_tab2[0]."_#".$contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
                $i++;
	}
	else
        {
		if ($contador>=2)
		{
			print $split_tab1[0]."_#".$contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		}
		elsif ($contador==1)		
		{
			print $split_tab1[0]."_#".$contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		}
		elsif ($i+1 > $#tab)
		{
			print $split_tab1[0]."_#1","\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		}
		
		#print $contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		#print $split_tab1[0]."_#".$contador,"\t",$split_tab1[1],"\t",$split_tab1[2],"\t",$split_tab1[3],"\n";
		$contador=1;
                $i++;
	}
}
