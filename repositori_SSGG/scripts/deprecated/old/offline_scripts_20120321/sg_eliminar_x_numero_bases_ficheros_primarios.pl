#!/usr/bin/perl -w

sub usage
{
        print "\nEste script elimina las 'x' últimas lecturas de un fichero csfasta (y su fichero .qual asociado).\n";
        print "\nCOMO SE USA: sg_eliminar_x_numero_bases_ficheros_primarios.pl <fichero_csfasta> <fichero_qual> <x>\n";
        print "ejemplo: sg_eliminar_x_numero_bases_ficheros_primarios.pl fichero_lecturas.csfasta fichero_calidades.qual 5 \n\n";
        print "INPUT: Dos ficheros y un entero: <fichero_csfasta> es el fichero con las lecturas, <fichero_qual> es el fichero con los valores de calidad asociados a las lecturas, y <x> es el número de lecturas a extraer de los ficheros de partida.\n\n";
        print "OUTPUT: Ficheros con la misma información que los de partida pero sustituyendo las 'x' últimas lecturas/valores por '.'. Los nombres de los ficheros son iguales a los de entrada pero añadiendo '_trimX' antes de su extensión.\n";
        exit(1);
}

# Si sólo ejecutamos el script, se imprime las instrucciones de uso
if(scalar(@ARGV) == 0)
{
        usage();
}

# Check arguments and extensions of the files
$extension_CSFASTA=substr($ARGV[0],length($ARGV[0])-8,8);
$extension_QUAL=substr($ARGV[1],length($ARGV[1])-5,5);

if (lc($extension_CSFASTA) ne ".csfasta" or lc($extension_QUAL) ne ".qual") {die "La extensión de los archivos de entrada -" . $extension_CSFASTA . ", " . $extension_QUAL . "- parece incorrecta!!";}

$name_CSFASTA=substr($ARGV[0],0,length($ARGV[0])-8);
$name_QUAL=substr($ARGV[1],0,length($ARGV[1])-5);

$x=int($ARGV[2]);

if ($x < 0) {die "El valor de lecturas a eliminar del fichero primario es incorrecto!!";}
print "Se procede a eliminar " . $x . " lecturas de los ficheros " . $name_CSFASTA . "\n";

# main body of the programm

# Abrimos los archivos de entrada
open(CSFASTA,"<",$ARGV[0]);
open(QUAL,"<",$ARGV[1]);

# the output files are opened in write-mode
open (OUTCSFASTA, ">", $name_CSFASTA . "_trimmed.csfasta");
open (OUTQUAL, ">", $name_QUAL . "_trimmed.qual");

# The 'csfasta' file will be parsed
while (my $line_csfasta = <CSFASTA>)
{
	chomp($line_csfasta);
	
	# Only the lines that not start with ">" are proccessed
	if (substr ($line_csfasta,0,1) ne ">")
	{
		if ($x < length($line_csfasta))
		{
			$line_csfasta=substr($line_csfasta,0,(length($line_csfasta)-$x));
		}
		#$line_csfasta=$line_csfasta . '.'x$x;
	}
		
	print OUTCSFASTA $line_csfasta,"\n";
}

# The 'qual' file will be parsed
while (my $line_qual = <QUAL>)
{
        chomp($line_qual);

        # Only the lines that not start with ">" are proccessed
        if (substr($line_qual,0,1) ne ">")
        {
		@split_line_qual=split(/ /,$line_qual);
		$NumElements=$#split_line_qual+1;
		
		$line_qual="";
	
		if (int($x) < $NumElements)
		{
			for ($i=0; $i<($NumElements-$x); $i += 1)
			{
				$line_qual=$line_qual . $split_line_qual[$i] . " ";
			}
		}
		chomp($line_qual);	
                #$line_qual=$line_qual . '.'x$x;
        }

        print OUTQUAL $line_qual,"\n";
}

