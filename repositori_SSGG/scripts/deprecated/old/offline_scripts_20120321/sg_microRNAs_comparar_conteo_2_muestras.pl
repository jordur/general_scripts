#!/usr/bin/perl
sub usage {
    print "\nPARA QUE SIRVE ESTO: Compara dos ficheros tabulados con el conteo de los miRNAs. Ej: OrR1 Vrs OrR2\n";
    print "Input: fichero .maduros de las dos muestras que se quieren comparar concatenados y ordenados por posicion de inicio\n";
    print "Output: hay que darle el output \n";
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
        if ((($valor_problema-$valor_inferior)<=8)  || ($valor_inferior == $valor_problema))
        {
                $contador=$contador+1;
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
		
		$contador=1;
                $i++;
	}
}
