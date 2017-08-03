#!/usr/bin/perl -w
# Selecciona la secuencias en may#####as de secuencias donde se mezclan fragmentos en may#####a y min#####a.
# El resultado lo imprime en un fichero fasta con el nombre secuencia.fasta

# Si no introducimos el ning#####chero, se imprime el modo de uso
sub usage {
        print "\n Selecciona solo las secuencias en mayusculas.\n";
        print "\nCOMO SE USA: sg_extraer_secuencia_en_mayusculas.pl <archivo>\n";
        print "Esta preparado para el formato\n\n";
        print "EXON 1 \n";
        print "aaatgctaACCTGGACTGAaaaaaa \n\n";
#       print "El resultado lo imprime en un fichero fasta con el nombre secuencia.fasta \n\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

# Definimos las variables
my $secuencia;
my @secuencia;
my $secuencia_final;

# Abrimos el archivo
open(SECUENCIA,"<",$ARGV[0]);
my $nombre;
# Lee la secuencia y no coge las l###as que empiecen por E
while  (my $secuencia_linea=<SECUENCIA>)
{
	$inicio = substr ($secuencia_linea,0,1);
	$secuencia_linea=~ s/ +/ /g;
	@split_linea= split (" ",$secuencia_linea);
	#my $nombre;
	if ($inicio eq ">")
	{
		$nombre = substr ($secuencia_linea,1,length($secuencia_linea)-2);
	}
	#print $split_linea[0],"\n";
	#$name=$nombre;
	else
	{	
		print $nombre,"_",$split_linea[0],"\t",$nombre,"\t",$split_linea[1],"\t",$split_linea[2],"\t",$split_linea[3],"\t",$split_linea[4],"\n";
	}


}
