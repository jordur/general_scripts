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
	$valor_final_inferior = $split_tab1[2];
	$valor_final_superior = $split_tab2[2];
        $valor_problema = $split_tab2[1];
	$valor_superior = $split_tab1[2];
	$valor_problema_superior = $split_tab2[2];
	if (($a<=$#tab) && ($valor_problema eq $valor_inferior) && ($valor_final_inferior eq $valor_final_superior))
        {
		#push (@number,$valor_inferior);
		#push (@number,$valor_problema);
		#push (@number,$valor_superior);
		#push (@number,$valor_problema_superior);
		
		$contador=$contador+1;
	}
	else
	{
		if ($contador >9)
		{
			#$max = $number[0];
    			#$min = $number[0];

    			#foreach $i (@number[1..$#number])
    			#{
        		#	if ($i > $max)
        		#	{
            		#		$max = $i;
        		#	}
        		#	elsif ($i < $min)
        		#	{
            		#		$min = $i;
        		#	}
    			#}
			print $split_tab1[0]."_x".$contador,"\t",$valor_inferior,"\t",$valor_final_inferior,"\t",$split_tab1[3],"\n";
		}
		$contador=1;
		#@number=();
	}
}
