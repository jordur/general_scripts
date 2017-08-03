#!/usr/bin/perl

sub usage
{
    print "\nCOMO SE USA:\n";
    print " Input: fichero de los miRNAs en dos columnas\n\t ";
    print " Output: fichero en formato fasta  \n\t ";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}

open(DOSCOLUMNAS,"<",@ARGV[0]);
open (FASTA, ">", 'blastclust.fasta') or die "No puedo abrir blastclust.fasta'.\n";

my @doscolumnas;
my @split_dos;
my $columnas_linea;
while  (my $columnas_linea=<DOSCOLUMNAS>)
{
        chomp($columnas_linea);
        push (@doscolumnas,$columnas_linea);
}


close (DOSCOLUMNAS);


for ($i =0;$i<=$#doscolumnas ; $i++)
{
	@split_dos = split (/\t/, $doscolumnas[$i]);
	print FASTA ">".$split_dos[0],"\n";
	print FASTA $split_dos[1],"\n";
}

close (FASTA);
