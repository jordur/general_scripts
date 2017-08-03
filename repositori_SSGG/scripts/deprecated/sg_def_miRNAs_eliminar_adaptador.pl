#!/usr/bin/perl

sub usage {
    print "\nPARA QUE SIRVE ESTO: Elimina los ultimos nt de la secuencia hasta obtener lecturas de 18nt\n";
    print "Input: fichero .csfasta con las lecturas no mapeadas\n";
    print "Output: fichero .csfasta con las lecturas mapeadas recortadas \n";
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
	$mayor = substr ($tab_linea,0,1);
	#print $mayor, "\n";
	if ($mayor ne '>')
	{
		$dieciocho_nt = substr ($tab_linea,0,18);
		print $dieciocho_nt, "\n";
	}
	else
	{
		print $tab_linea,"\n";
	}
}

print $contador,"\n";
close (TAB);

