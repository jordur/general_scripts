#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: coloca en columnas los valores de cuantificacion para cada muestra en estudio. Si el valor no existe se sustituye por FALTA. \n";
    print " Input: 1) fichero tabulado del conteo de la muestra 2) fichero con el nombre de los genes que existen y su longitud  3) numero de columna del conteo (se asume que la primera columna corresponde siempre al nombre del gen) 4) numero de columna con la longitud transcrito \n\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(CONTEOPORGEN,"<",@ARGV[0]);
open(GENEUNICOS,"<",@ARGV[1]);
my $numero_columna=@ARGV[2];
my $longitud_transcrito=@ARGV[3];
#open(SENSE,">","$nombre_chr.conocidos.sense.ma.tab") or die "No puedo abrir $nombre_chr.conocidos.sense.ma.tab.\n";
#open(ANTISENSE,">","$nombre_chr.conocidos.antisense.ma.tab") or die "No puedo abrir $nombre_chr.conocidos.antisense.ma.tab.\n";
#open(NUEVOS,">","$nombre_chr.nuevos.ma.tab") or die "No puedo abrir $nombre_chr.nuevos.ma.tab.\n";


#my @conteo;
#my @mirbase2;

while  (my $conteo=<CONTEOPORGEN>)
{
        chomp($conteo);
	$almohadilla=substr($conteo,0,1);
	if($almohadilla ne "#")
	{
        	push (@conteoporgen,$conteo);
	}
}

close (CONTEOPORGEN);

while  (my $unico=<GENEUNICOS>)
{
        chomp($unico);
        push (@geneunicos,$unico);
}

close (@geneunicos);
#print join ("\n",@conteoporgen),"\n";
#print $#geneunicos;
$contador=0;

for (my $a=0; $a<=$#geneunicos ; $a++)
{
	@split_geneunicos = split ("\t",$geneunicos[$a]);
	for (my $i=0; $i<=$#conteoporgen;$i++)
	{
		@split_conteoporgen= split ("\t",$conteoporgen[$i]);
		if ($split_conteoporgen[0] eq $split_geneunicos[0])
		{
			print $geneunicos[$a],"\t",$split_conteoporgen[$longitud_transcrito-1],"\t",$split_conteoporgen[$numero_columna-1],"\n";
			$contador=1;
			last;
		}
	}
	if ($contador==0)
	{
		print $geneunicos[$a],"\tFALTA\n";
	}
	$contador=0;
}
