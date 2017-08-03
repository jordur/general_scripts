#!/usr/bin/perl -w
# Selecciona la secuencias en mayúsculas de secuencias donde se mezclan fragmentos en mayúscula y minúscula.
# El resultado lo imprime en un fichero fasta con el nombre secuencia.fasta

# Si no introducimos el ningún fichero, se imprime el modo de uso
sub usage {
    	print "\n Selecciona solo las secuencias en mayusculas.\n";
	print "\nCOMO SE USA: sg_extraer_secuencia_en_mayusculas.pl <archivo>\n";
    	print "Esta preparado para el formato\n\n";
    	print "EXON 1 \n";
    	print "aaatgctaACCTGGACTGAaaaaaa \n\n";
#	print "El resultado lo imprime en un fichero fasta con el nombre secuencia.fasta \n\n";
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

# Lee la secuencia y no coge las líneas que empiecen por E
while  (my $secuencia_linea=<SECUENCIA>)
{	
	my $primer_caracter = substr ($secuencia_linea,0,1);
#	print  $primer_caracter,"\n";
	if ($primer_caracter ne "E")
        {
		chomp($secuencia_linea);
        	#print $secuencia_linea,"\n";
		$secuencia_linea =~ s/a//g;
		$secuencia_linea =~ s/t//g;
		$secuencia_linea =~ s/g//g;
		$secuencia_linea =~ s/c//g;
		$secuencia_linea =~ s/n//g;
		#print $secuencia_linea,"\n";		
		# Eliminamos todos los espacios
		$secuencia_linea =~ s/\s//g;
        	push (@secuencia,$secuencia_linea);
	}
}
close (SECUENCIA);

#Unimos todas las secuencias en un línea
print join ('',@secuencia),"\n";

# Imprimimos los resultados en un fichero fasta llamado secuencia
#$secuencia_final = "secuencia.fasta";
#open (OUTPUT,">$secuencia_final");
#print OUTPUT ">secuencia";
#print OUTPUT $secuencia_final;
#close (OUTPUT);

#Salimos del programa
exit;

