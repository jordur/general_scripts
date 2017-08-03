#! /usr/bin/perl -w

# Descripción
# Script diseñado para sustituir una base heterocigota en la base homocigota diferente a la referencia. Está hecho con una subrutina, por lo que se puede pegar a cualquier script

#use strict;

sub usage
{
        print "\nCOMO SE USA: sg_convertir_SNPs_de_heterocigotos_a_homocigotos.pl >input>\n\n";
        print "INPUT: chr1    45162   rs10399749      C       T,A     0       2\n\n";
        print "OUTPUT: chr1    45162   rs10399749      C       M \n\n";
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
	my @alternative = split (",",$ex[4]);
	if (($ex[5] eq "0") && ($ex[6] eq "1"))
	{
		$base=$ex[3].$alternative[0];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";	
	}
	elsif (($ex[5] eq "0") && ($ex[6] eq "2"))
	{
		$base=$ex[3].$alternative[1];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "1") && ($ex[6] eq "0"))
	{
		$base=$alternative[0].$ex[3];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "1") && ($ex[6] eq "1"))
	{
		$base=$alternative[0];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "1") && ($ex[6] eq "2"))
	{
		$base=$alternative[0].$alternative[1];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "2") && ($ex[6] eq "0"))
	{
		$base=$alternative[1].$ex[3];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "2") && ($ex[6] eq "1"))
	{
		$base=$alternative[1].$alternative[0];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "2") && ($ex[6] eq "2"))
	{
		$base=$alternative[1];
		$a= convertir_a_homocigoto ($base);
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
	elsif (($ex[5] eq "0") && ($ex[6] eq "0"))
	{
		$a=$ex[3];
		print $ex[0],"\t",$ex[1],"\t",$ex[2],"\t",$ex[3],"\t",$a,"\n";
	}
 
}

close(SNP);


sub convertir_a_homocigoto 
{
        if(($base eq "AG") || ($base eq "GA"))
	{
		return "R";
	}
	elsif(($base eq "AC") || ($base eq "CA") )
	{
                return "M";
	}
        elsif(($base eq "AT") || ($base eq "TA") )
        {
                return "W";
	}
	elsif(($base eq "GC") || ($base eq "CG") )
        {
                return "S";
        }
	elsif(($base eq "CT") || ($base eq "TC"))
        {
                return "Y";
        }
	elsif(($base eq "GT") || ($base eq "TG"))
        {
                return "K";
        }
	else
	{
		return $base;
	}
	
}

