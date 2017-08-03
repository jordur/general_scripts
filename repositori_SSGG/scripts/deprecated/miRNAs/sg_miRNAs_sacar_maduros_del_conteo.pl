#!/usr/bin/perl
sub usage
{
    print "\nPARA QUE SIRVE: compara fichero del conteo de miRBase con el fichero de miRNAs maduros\n";
    print " Input: 1) fichero tabulado del conteo de miRNAs 2) fichero con miRNAs maduros ej: mir2<TAB>chr<TAB>inicio<TAB>final<TAB>orientacion HAY QUE DARLE EL NOMBRE DE LA SALIDA >nombre_fichero_salida\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(CONTEO,"<",@ARGV[0]);
open(MIRBASE,"<",@ARGV[1]);

my @conteo;
my @mirbase2;

while  (my $linea2=<CONTEO>)
{
        chomp($linea2);
        push (@conteo,$linea2);
}

close (CONTEO);


while  (my $linea=<MIRBASE>)
{
        chomp($linea);
        push (@mirbase2, $linea);
}

close (MIRBASE);

my $contador = 0;
for (my $a=0; $a<=$#mirbase2 ;$a++)
{
	@mirbase= split ("\t",$mirbase2[$a]);	
	for (my $i=0;$i<=$#conteo;$i++)
	{
		@conteo_split= split ("\t",$conteo[$i]);
		@conteo_split2 = split ("_x",$conteo_split[0]);
		if (((abs($conteo_split[1]-$mirbase[2])) <= 10))
		{
			print $mirbase2[$a],"\t",$conteo_split2[1],"\n";
		}
	}
}
