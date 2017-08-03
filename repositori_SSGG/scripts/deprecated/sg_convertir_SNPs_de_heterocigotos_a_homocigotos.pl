#! /usr/bin/perl -w

# Descripción
# Script diseñado para sustituir una base heterocigota en la base homocigota diferente a la referencia. Está hecho con una subrutina, por lo que se puede pegar a cualquier script

use strict;

sub usage
{
        print "\nCOMO SE USA: sg_convertir_SNPs_de_heterocigotos_a_homocigotos.pl >input>\n\n";
        print "EJEMPLO: sg_convertir_SNPs_de_heterocigotos_a_homocigotos.pl SNPs_no_descritos.txt\n\n";
        print "INPUT: El fichero de input debe ser un fichero tabulado con las siguientes columnas:\ncolumna 1: nombre del cromosoma (Sólo se admite el formato 'chr1')\ncolumna 2: la posición en el cromosoma.\ncolumna 3: base en la referencia.\ncolumna 4: base detectada.\ncolumna 5: calidad del consenso.\ncolumna 6: calidad del SNP.\ncolumna 7: calidad máxima de mapeo.\ncolumna 8: coverage. Si el fichero de input no tiene de las columnas 5 a 8, el script funciona igual.\n\n";
        print "OUTPUT: Obtenemos un archivo igual al input donde se han sustituído las bases heterocigotas por la base homocigota diferente a la referencia. No acepta las bases IUPAC V, H, D, B, X o N.\n\n";
        exit(1);
}

if(scalar(@ARGV) == 0){
    usage();
}

# Abrimos el archivo de SNPs

open(SNP,'<',$ARGV[0]);

while (my $snp=<SNP>)
{
	chomp($snp);
	my @ex = split(/\t/,$snp);
	my $base = convertir_a_homocigoto($ex[2],$ex[3]);
	print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$base,"\t",$ex[4],"\t",$ex[5],"\t",$ex[6],"\t",$ex[7],"\n";
}

close(SNP);

# Salimos del programa
exit;

##########################################################################################
#
# Subrutina convertir a homocigoto
#
##########################################################################################

# Necesita dos variables. La primera debe ser la referencia y la segunda el genotipo. No acepta IUPAC V, H, D, B, X o N

sub convertir_a_homocigoto
{
	my $out;
        if(($_[1] eq "A") || ($_[1] eq "C") || ($_[1] eq "G") || ($_[1] eq "T"))
        {
		return $_[1];
        }
        else
        {
                if($_[1] eq "R")
		{
			if($_[0] eq "A")
			{
				$out = "G";
			}
			else
			{
				$out = "A";
			}
			return $out;
		}
		elsif($_[1] eq "M")
		{
			if($_[0] eq "A")
                        {
                                $out = "C";
                        }
                        else
                        {
                                $out = "A";
                        }
                        return $out;
		}
                elsif($_[1] eq "W")
                {
                        if($_[0] eq "A")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "A";
                        }
                        return $out;
		}
		elsif($_[1] eq "S")
                {
                        if($_[0] eq "C")
                        {
                                $out = "G";
                        }
                        else
                        {
                                $out = "C";
                        }
                        return $out;
                }
		elsif($_[1] eq "Y")
                {
                        if($_[0] eq "C")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "C";
                        }
                        return $out;
                }
		elsif($_[1] eq "K")
                {
                        if($_[0] eq "G")
                        {
                                $out = "T";
                        }
                        else
                        {
                                $out = "G";
                        }
                        return $out;
                }
		else
		{
		die("este script no acepta IUPAC V, H, D, B, X o N");
		}
	}
}

