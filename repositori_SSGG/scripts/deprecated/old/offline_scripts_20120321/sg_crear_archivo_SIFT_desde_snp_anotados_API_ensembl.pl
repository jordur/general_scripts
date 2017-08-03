#!/usr/bin/perl -w

sub usage {
	print "\nEste script crea un archivo en formato SIFT del fichero tabulado que obtenemos de la anotación con la API del ensembl\n";
        print "\nCOMO SE USA: sg_crear_archivo_SIFT.pl <fichero> \n";
        print "ejemplo: sg_crear_archivo_SIFT.pl snp_no_descritos_chr2.txt\n\n";
	print "INPUT: fichero tabulado en el que la primera coordenada es el cromosoma, la segunda la coordenada, la tercera, la base referencia, la cuarta, el genotipo encontrado y la doce, la cadena.\n\n";
	print "OUTPUT: fichero con el formato 2,3687911,1,T/G y cuyo nombre es igual al del input con el prefijo SIFT_\n\n";
        exit(1);

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
}
if(scalar(@ARGV) == 0){
    usage();
}

# Abrimos el archivo de salida
open(OUTPUT,">","SIFT");

# Abrimos el archivo 
open(FILE,"<",$ARGV[0]);

while (my $snv=<FILE>)
{
        chomp($snv);
	my @ex = split(/\t/,$snv);
	if($ex[1] ne "coordinate")
	{
		$ex[0] =~ tr/chr//d;
		$genotype = convertir_a_homocigoto($ex[2],$ex[3]);
		print OUTPUT "$ex[0],$ex[1],$ex[11],$ex[2]/$genotype\n";
	} 
}

close(FILE);

# Cerramos los ficheros de salida
close(OUTPUT);

system("sort -u SIFT > SIFT_$ARGV[0]");
system("rm SIFT");

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

##############################################################################3
