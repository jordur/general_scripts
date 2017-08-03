#!/usr/bin/perl
# Elimina el primer cero y el espacio siguiente de los ficheros de calidad resultante de convertir los archivos fastq con el script fastq2csfasta.pl


sub usage
{
    print "\nPARA QUE SIRVE:  Elimina el primer cero y el espacio siguiente de los ficheros de calidad resultante de convertir los archivos fastq con el script fastq2csfasta.pl \n COMO SE USA:\n";
    print "sg_eliminar_primer_cero_qual_fastq.pl <input> > <output>\n";
    print "Input: fichero de calidad resultante de la conversión del archivo fastq.\n";
    print "Output: Hay que darle un fichero de salida. \n";
    print "\n";
    exit(1);
}

if(scalar(@ARGV) == 0)
{
    usage();
}


open(QUAL,"<",@ARGV[0]);

while (my $linea = <QUAL>)
{
	chomp ($linea);
	$mayor= substr($linea,0,1);
	if ($mayor eq ">")
	{
		print $linea,"\n";	
	}
	else
	{
		$qv = substr($linea,2, length $linea);
		print $qv,"\n";
	}
}

close (QUAL);
