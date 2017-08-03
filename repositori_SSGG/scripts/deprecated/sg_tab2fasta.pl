#!/usr/bin/perl

sub usage
{
    print "\nPARA QUE SIRVE: transforma un fichero tabulado en un fichero fasta \nCOMO SE USA:\n";
    print " Input: 1) fichero tabulado 2) columna que contiene el nombre de la secuencia 3) columna que contiene la secuencia o los valores de calidad \n";
    print " Output: capturar porque sale en pantalla  \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


open(FICHERO,"<",@ARGV[0]);
$col_nombre = @ARGV[1];
$col_secuencia = @ARGV[2];

while (my $linea = <FICHERO>)
{
	chomp ($linea);
	@split_linea = split ("\t",$linea);
	print ">",$split_linea[$col_nombre-1],"\n";
	print $split_linea[$col_secuencia-1],"\n";
}

close (FICHERO);
