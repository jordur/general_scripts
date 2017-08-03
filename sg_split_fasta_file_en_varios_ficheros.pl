#!/usr/bin/perl
sub usage {
    print "\nPARA QUE SIRVE: Divide un fichero en formato FASTA varios ficheros fasta con un numero de secuencias determinado\n";
    print "COMO SE USA:\n";
    print "sg_split_fasta_file_en_varios_ficheros.pl fichero_fasta numero_secuencias_por_fichero\n";
    print "\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

open(FASTA,"<",@ARGV[0]);
my $desc;
$num = @ARGV[1]-1;
$contador=0;
$con=0;
my $lineas;
my $desc;
chomp($lineas);
while  ($lineas=<FASTA>)
{
        chomp($lineas);
	$primerChar = substr ($lineas,0,1);
	if ($contador <= $num)
	{
		push (@array, $lineas);
		if ($primerChar eq ">")
		{
		    if ($lineas =~ />lcl\|(.+)/){
			$desc = $1;
			#print $desc,"\n";
		    }
		    #&escribir_fichero(@array),"\n";
			$contador = $contador+1;
			
		}
	}
	elsif (($contador > $num))
	{
		
		if ($primerChar eq ">")
		{
		    if ($lineas =~ />lcl\|(.+)/){
			$desc= $1;
		    }
		    &escribir_fichero(@array),"\n";
		    @array=();
		    push (@array, $lineas);
		    $contador=1;
		}
	else{
		    push (@array, $lineas);
	    }
	}
	
	
}
#Imprimo de nuevo el array porque si llega al final del fichero y el contador no es igual que el numero de secuencias que se ha puesto como input las ultimas lineas no se escriben.
&escribir_fichero(@array),"\n";

sub escribir_fichero
{
	@lalinea = @_;
	$con = $con+1;
	#print $desc,"\n";
	open (OUT,">","resultado_".$con);
	#open (OUT, ">",$desc."_".$con);
	print OUT join ("\n",@lalinea),"\n";
	close (OUT)
}
close (FASTA);
