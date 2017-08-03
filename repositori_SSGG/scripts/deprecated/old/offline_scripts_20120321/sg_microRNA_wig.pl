#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE VALE:\n";
    print "Input: 1) Fichero tabulado del conteo de cada grupo de miRNAs y 2) Cromosoma sin cabecera  3) Nombre del cromosoma ej:2L,3R\n\t ";
    print "Output: Fichero wig \n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(GRUPOS,"<",@ARGV[0]);

while (my $linea = <GRUPOS>)
{
	chomp ($linea);
	@array0=split ("\t",$linea);
	@array1=split ("x",$array0[0]);
	push (@uno,$array0[1]);
	push (@dos,($array0[2]-$array0[1]));
	push (@conteo,$array1[1]);
	
}

open(REFERENCIA,"<",@ARGV[1]);
while  (my $lineas=<REFERENCIA>)
{
        chomp($lineas);
        $secuencia.=$lineas;
}

close (GRUPOS);
close (REFERENCIA);

$longitud = length ($secuencia);
@array=(1..$longitud);
$chr=$ARGV[2];
print "browser position chr$chr:$uno[0]-",$uno[-1]+$dos[-1],"\nbrowser hide all\ntrack type=wiggle_0 name=\"SOLiD genome coverage\" description=\"SOLiD genome coverage\" visibility=full color=255,0,0 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=chr$chr span=1\n";

for ($i=0; $i<=$longitud ; $i++)
{
	push (@posicion,"0");
}	

for ($a=0;$a<=$#uno; $a++)
{
	for ($meter=$uno[$a]; $meter<=$uno[$a]+$dos[$a]; $meter++)
	{
		$nuevo_valor=$posicion[$meter-1]+$conteo[$a];
		splice (@posicion,$meter-1,1,$nuevo_valor);
	}
	
}
for ($e=0;$e<=$#posicion;$e++)
{
	if ($posicion[$e]>0)
	{
		print $array[$e],"\t",$posicion[$e],"\n";
	}
}
