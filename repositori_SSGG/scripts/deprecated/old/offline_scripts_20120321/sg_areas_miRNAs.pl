#!/usr/bin/perl
use List::Util qw (sum);

my @fila;#Del coverage.txt.
my @fila2;
my @lineas;

#print "\nCoverage minimo: ";
#my $coverage_maximo= <STDIN>;
#print "\nNumero minimo de bases entre areas: ";
#my $numero_bases_entre_areas = <STDIN>;
#print "\nCalculando high coverage ...\n";


#Suma al array SUMA.

#open (FILE1, "<", 'coverage.txt');
#open (FILE1, "<", @ARGV);#Fichero de entrada como parametro. I/O en ruta local usuario.
$string1= '34	4';
$string2= '35	10';
$string3= '40	3';
$string4= '41	5';
push (@lineas, $string1,$string2, $string3, $string4);

for (my $i=0; $i<=$#lineas+1; $i++)
{
	
	@fila = split (/\s+/,@lineas[$i]);
	@fila = @fila2;	

}

print @fila[0],"\n";


#print @fila[1];
#close (FILE1);#Cierro la entrada del fichero coverage.txt


